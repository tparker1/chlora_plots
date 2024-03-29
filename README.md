# WIP 
Please visit again later

## Step 1: masking/create_masks.ipynb

### Initialize
Initialize the masking process with an encompassing bounding box (of min/max lon/lat indices). 
_Give **indices** for a box that is at least larger than all of the polygons you are interested in._

Example: for the "Arrigo et al 2017 Sample Areas.shp" regions, I use initial bounds:

minx_index, maxx_index, miny_index, maxy_index = "11800", "13029", "2894", "2325"

The first_day, last_day is arbitrary. 

_Note: You will need to figure out which indices correspond to the initial encompassing bounding box._

### load polygons
SET polygon_folder  (path containing polygon_file)
SET polygon_file    (a single shp file with one to many polygons, such as "Arrigo et al 2017 Sample Areas.shp")

SET buffer_distance to expand the polygons

#### Subset polygons
After running the cell to load in all polygons in polygon_file, create subset (group) if needed. 

SET regions - a list of Region IDs. 

Subsetting is not necessary, but setting group as a group identifer is.  

### Create Masks
_This will store masked arrays on your system as pkl files (one .pkl for each polygon region). You should CHANGE the output path to an appropriate folder on your system._
   
The big idea:

The given polygons (ex: Arrigo et al 2017 Sample Areas.shp) are mapped to the Ocean Colour lat/lon grid and saved (via a pkl file) as an numpy masked array.

_IMPORTANT NOTE: The code is set up to run on ALL polygons in the shp file. If you only want a subset, change the "for x in range(len(polygons))" to something like "for x in x_of_interest" where "x_of_interest" is a list of integers corresponding to the polygon indicies you are interested in._

_Example: To only get masks for Region 1 and Region 4 in "Arrigo et al 2017 Sample Areas.shp", I would set "x_of_interest = [1,3]"._

### Plot
Can be used to visually check that your polygon masks are where you expect.

## Step 2: download/retrieve_opendap_gl_chla.ipynb

### Load Regional Masks
In step 1 you created and stored numpy masks for each region as a .pkl file. Now we will load them back into memory. CHANGE path_to_pkls_folder such that it points to the folder containing your pkl files. 

A
_Note: If your masked files are not sequentially named then update the for x in range(6)_






