import os
import numpy as np
import netCDF4 as nc

import matplotlib.pyplot as plt
import matplotlib.patheffects as pe

import geopandas as gpd
import cartopy.crs as ccrs
from shapely.geometry import Point, shape, Polygon

from scipy.ndimage import center_of_mass




def get_gl_extents():
    # return {'min_lat': 58.7447917, 'max_lat': 88.953125, 'min_lon': -70.869792, 'max_lon': -8.1197917}
    # return {'min_lat': 58, 'max_lat': 86, 'min_lon': -80, 'max_lon': -7}
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



def expand_bounds(polygon, distance):
    return polygon.buffer(distance)


def load_polygons(polygon_filepath, buffer_distance=0):

    polygons  = gpd.GeoDataFrame.from_file(polygon_filepath)

    polygons['geometry'] = polygons['geometry'].apply(expand_bounds, distance=buffer_distance)
    return polygons



def save_polygon_group(group, polygon_folder, filename, regions, polygons):
    group_folder = os.path.join(polygon_folder, f'group_{group}')

    # Write the polygon group to a shapefile
    if not os.path.exists(os.path.join(polygon_folder, group_folder)):
        os.makedirs(os.path.join(polygon_folder, group_folder))
    polygon_filepath = os.path.join(polygon_folder, group_folder, filename+'.shp')

    # Save a text file with the region list
    region_file = os.path.join(polygon_folder, group_folder, filename+'.txt')
    with open(region_file, 'w') as f:
        for region in regions:
            f.write(region + '\n')

    polygons.to_file(os.path.join(polygon_filepath), driver='ESRI Shapefile')

    print(f"Saved {filename} to {polygon_filepath}")


def create_gl_polygon():
    gl_extents = get_gl_extents()

    gl_polygon = Polygon([(gl_extents['min_lon'], gl_extents['min_lat']),
                        (gl_extents['max_lon'], gl_extents['min_lat']),
                        (gl_extents['max_lon'], gl_extents['max_lat']),
                        (gl_extents['min_lon'], gl_extents['max_lat'])])

    # Create a GeoDataFrame with the polygon geometry
    return gpd.GeoDataFrame(geometry=[gl_polygon])



def visualize_bounding_box(lons, lats, bounding_box):
    gl_extents = get_gl_extents()
    
    minx_index, maxx_index, miny_index, maxy_index = bounding_box
    
    fig, ax = plt.subplots(figsize=(10, 10), dpi=100, subplot_kw={'projection': ccrs.Stereographic(central_longitude=-45, central_latitude=90)})

    ax.set_extent([gl_extents['min_lon']+10, gl_extents['max_lon']-10, gl_extents['min_lat'], gl_extents['max_lat']-1], crs=ccrs.PlateCarree())

    # plot the bounding box as a shaded rectangle
    ax.fill([lons.min(), lons.max(), lons.max(), lons.min(), lons.min()],
            [lats.min(), lats.min(), lats.max(), lats.max(), lats.min()],
            transform=ccrs.PlateCarree(), color='green', alpha=0.5)

    ax.coastlines()

    # Print some reference info about your box
    title = "Bounding Box Information"
    lat_lon_info = "Min Lon: {:<10.5f} | Max Lon: {:<10.5f} | Min Lat: {:<10.5f} | Max Lat: {:<10.5f}".format(lons.min(), lons.max(), lats.min(), lats.max())
    index_info = "MinX Index: {:<6} | MaxX Index: {:<6} | MinY Index: {:<6} | MaxY Index: {:<6}".format(minx_index, maxx_index, miny_index, maxy_index)

    title_fontsize = 18
    text_placement = +0.01 # Use this to adjust the vertical placement of the text

    plt.figtext(0.5, fig.subplotpars.top + text_placement, title, ha='center', fontsize=title_fontsize)
    plt.figtext(0.5, fig.subplotpars.bottom + text_placement-0.04, lat_lon_info, ha='center', fontsize=title_fontsize-6)
    plt.figtext(0.5, fig.subplotpars.bottom + text_placement-.08, index_info, ha='center', fontsize=title_fontsize-6)

    plt.show()
    plt.close()




def plot_polygons(polygons, lons, lats):
    gl_extents = get_gl_extents()

    fig, ax = plt.subplots(figsize=(10, 10), dpi=100, subplot_kw={'projection': ccrs.Stereographic(central_longitude=-45, central_latitude=90)})

    ax.set_extent([gl_extents['min_lon']+20, gl_extents['max_lon']-20, gl_extents['min_lat'], gl_extents['max_lat']-1], crs=ccrs.PlateCarree())

    # plot the bounding box as a shaded rectangle
    ax.fill([lons.min(), lons.max(), lons.max(), lons.min(), lons.min()],
            [lats.min(), lats.min(), lats.max(), lats.max(), lats.min()],
            transform=ccrs.PlateCarree(), color='green', alpha=0.2)
    
    for index, row in polygons.iterrows():
        geometry = row['geometry']
        plons, plats = geometry.exterior.xy
        ax.fill(plons, plats, transform=ccrs.PlateCarree(), color='green', alpha=0.5)
        ax.plot(plons, plats, transform=ccrs.PlateCarree(), color='green', alpha=0.8, linewidth=.5)


    ax.set_title('Bounding Box \nwith Buffered Regions', fontsize=14)
    ax.coastlines()

    plt.show()



def visualize_masks(masked_arrays, lons, lats, path_to_pkls_folder, group):
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

    ax.set_title('Region Masks', fontdict={'fontsize': 20})
    ax.coastlines()

    # Save the plot to the mask folder
    mask_plot = os.path.join(path_to_pkls_folder, f'mask_plot_group_{group}.png')
    fig.savefig(mask_plot)

    plt.show()
    plt.close()







