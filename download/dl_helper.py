import os
import re
import shutil

import time as tm

import numpy as np

import pandas as pd
import netCDF4 as nc

import matplotlib.pyplot as plt
import matplotlib.patheffects as pe

import cartopy.crs as ccrs

from joblib import Parallel, delayed


def get_gl_extents():
    return {'min_lat': 58, 'max_lat': 89, 'min_lon': -81, 'max_lon': -3}

def get_lons_lats(bounding_box):
    minx_index, maxx_index, miny_index, maxy_index = bounding_box

    firstday, lastday = "4348", "4350"

    url = "https://www.oceancolour.org/thredds/dodsC/CCI_ALL-v6.0-1km-DAILY?lat["+maxy_index+":1:"+miny_index+"],lon["+minx_index+":1:"+maxx_index+"],chlor_a["+firstday+":1:"+lastday+"]["+maxy_index+":1:"+miny_index+"]["+minx_index+":1:"+maxx_index+"],time["+firstday+":1:"+lastday+"]"
    ds = nc.Dataset(url)

    lons = ds.variables['lon'][:]
    lats = ds.variables['lat'][:]

    ds.close()
    return lons, lats


def get_url(i, step, bounding_box):
    minx_index, maxx_index, miny_index, maxy_index = bounding_box 
    firstday = str(i)
    lastday = str(i+(step-1))
    url = "https://www.oceancolour.org/thredds/dodsC/CCI_ALL-v6.0-1km-DAILY?lat["+maxy_index+":1:"+miny_index+"],lon["+minx_index+":1:"+maxx_index+"],chlor_a["+firstday+":1:"+lastday+"]["+maxy_index+":1:"+miny_index+"]["+minx_index+":1:"+maxx_index+"],time["+firstday+":1:"+lastday+"]"
    return nc.Dataset(url)



def calculate_mean(mask_name, mask, t, step, time, bounding_box, max_chl=50):
    data = []

    ds = get_url(t, step, bounding_box)
    for t in range(len(time)):

        #chaotic progress bar
        if mask_name == '007':
            if (t+1) % max(1, len(time) // 20) == 0:  # Check if the current iteration is a multiple of 5% of the total number of iterations
                print(f"Step Progress: {t+1:6f} of {len(time):6f}", end="\r")

        chlor_a = np.ma.masked_where(mask.mask, ds.variables['chlor_a'][t].data)    
        chlor_a = np.ma.masked_outside(chlor_a, 0, max_chl)

        mean = chlor_a.mean()
        data.append({'region': mask_name, 'time': time[t], 'mean': mean})

    return pd.DataFrame(data)



def process_mask_index(start, step, masked_arrays, bounding_box, max_chla, destination_folder):

    for t in range(start, 10000, step):
        
        progress = t/10000 * 100
        print(f"Overall progress: {progress:3.0f}ish% complete")

        ds = get_url(t, step, bounding_box)    
        time = ds.variables['time'][:]

        results = Parallel(n_jobs=-1, backend="threading")(delayed(calculate_mean)(mask_name, masked_array, t, step, time, bounding_box, max_chla) for mask_name, masked_array in masked_arrays.items())

        df = pd.concat(results)
        
        output_file = os.path.join(destination_folder,  f'mean_chl_{t}_{t+(step-1)}.csv')
        df.to_csv(output_file, index=False)
        



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                            # magic
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# innards of 
# retrieve_opendap_gl_chla_data
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def get_ds(firstday, lastday, bounding_box):
    minx_index, maxx_index, miny_index, maxy_index = bounding_box

    url = "https://www.oceancolour.org/thredds/dodsC/CCI_ALL-v6.0-1km-DAILY?lat["+maxy_index+":1:"+miny_index+"],lon["+minx_index+":1:"+maxx_index+"],chlor_a["+firstday+":1:"+lastday+"]["+maxy_index+":1:"+miny_index+"]["+minx_index+":1:"+maxx_index+"],time["+firstday+":1:"+lastday+"]"
    ds = nc.Dataset(url)

    return ds

def get_mmddyy(timesinceepoch):
    return tm.gmtime(timesinceepoch * 86400)

def process_batch(i, step, bounding_box, lons, lats, destination_path, staging_path):
    firstday = str(i)
    lastday = str(i + step - 1)
    ds = get_ds(firstday, lastday, bounding_box)
    save_batch(ds, firstday, lons, lats, destination_path, staging_path)
    ds.close()

def save_batch(ds, firstday, lons, lats, destination_path, staging_path):
    time_var = ds.variables['time'][:]
    
    year = get_mmddyy(time_var[0]).tm_year
    print(f"Working on year {year}...", end="\r")

    fp = os.path.join(staging_path, f'chlor_a_data_{firstday}.nc')
    with nc.Dataset(fp, 'w') as f:
        f.createDimension('time', len(time_var))
        f.createDimension('lat', len(lats))
        f.createDimension('lon', len(lons))

        time_var_out = f.createVariable('time', 'f8', ('time',))
        lat_out = f.createVariable('lat', 'f8', ('lat',))
        lon_out = f.createVariable('lon', 'f8', ('lon',))
        chlor_a_out = f.createVariable('chlor_a', 'f8', ('time', 'lat', 'lon'))

        time_var_out[:] = time_var
        chlor_a_out[:] = ds.variables['chlor_a'][:]
        lat_out[:] = lats
        lon_out[:] = lons

    # move the file to the destination folder
    shutil.move(fp, os.path.join(destination_path, f'chlor_a_data_{firstday}.nc'))

def get_steps_remaining(destination_path):
    start, stop, step = 0, 9600, 50

    destination_files = [f for f in os.listdir(destination_path) if f.endswith('.nc')]

    completed_steps = [int(f.split('.')[0].split('_')[-1]) for f in destination_files]

    steps_needed = [i for i in range(start, stop, step) if i not in completed_steps]

    return step, steps_needed
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

def get_latest_file_number(destination_path):
    files = os.listdir(destination_path)
    files = [f for f in files if f.endswith('.csv')]

    file_numbers = [int(re.search(r'(\d+)\D*$', f).group(1)) for f in files]
    
    return max(file_numbers, default=-1) + 1



def visualize_masks(masked_arrays, lons, lats,): #, bounding_box):
    gl_extents = get_gl_extents()

    # Create a map in the EPSG 3413 projection
    fig, ax = plt.subplots(figsize=(10, 10), dpi=100, subplot_kw={'projection': ccrs.Stereographic(central_longitude=-45, central_latitude=90)})

    # Set the limits of the x and y axes
    ax.set_extent([gl_extents['min_lon']+20, gl_extents['max_lon']-20, gl_extents['min_lat'], gl_extents['max_lat']-1], crs=ccrs.PlateCarree())
    
    for mask_name, masked_array in masked_arrays.items():

        # Overlay the masked array on the map
        ax.imshow(masked_array, origin='upper', extent=[lons.min(), lons.max(), lats.min(), lats.max()], transform=ccrs.PlateCarree(), cmap='Accent', alpha=0.8)
       
        # Calculate the center of the unmasked portion (for fancy labeling)
        unmasked_indices = np.where(~np.isnan(masked_array))
        center_lat = lats[unmasked_indices[0]].mean()
        center_lon = lons[unmasked_indices[1]].mean()

        # ATTENTION change ha or va depending on where your regions are 
        ax.text(center_lon, center_lat, u"\u2190" + mask_name, ha='left', va='center', fontsize=10, fontweight='bold', color='black',
            path_effects=[pe.withStroke(linewidth=3, foreground='white')], transform=ccrs.PlateCarree())

   # plot the bounding box as a shaded rectangle
    ax.fill([lons.min(), lons.max(), lons.max(), lons.min(), lons.min()],
            [lats.min(), lats.min(), lats.max(), lats.max(), lats.min()],
            transform=ccrs.PlateCarree(), color='green', alpha=0.2)

    ax.set_title('Region Masks and Bounding Box', fontdict={'fontsize': 20})
    ax.coastlines()

    plt.show()
    plt.close()