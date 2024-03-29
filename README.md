# WIP 
Please visit again later

## Step 1: masking/create_masks.ipynb
[masking/create_masks.ipynb](https://github.com/tparker1/chlora_plots/blob/c0e34ab347d5392280ceb00c8096bde6909150a6/masking/create_masks.ipynb) 

### Initialize
Initialize the masking process with an encompassing bounding box (of min/max lon/lat **indices**). 
_Give **indices** for a box that is at least larger than all of the polygons you are interested in._

Example: for the _Arrigo et al 2017 Sample Areas.shp_ regions, I use initial bounds:

>minx_index, maxx_index, miny_index, maxy_index = "11800", "13029", "2894", "2325"

``visual check`` Plots the current bounding box on a map of Greenland.


### load polygons
``SET polygon_folder``  (path containing polygon_file)

``SET polygon_file``    (a single shp file with one to many polygons, such as _Arrigo et al 2017 Sample Areas.shp_)

``SET group`` name **_note: if you subset the polygons, you must set group name again_**

``SET buffer_distance`` to expand the polygons

#### Subset polygons
Not necessary, but useful if you do not want to run all polygons at once. 
``SET regions`` - a list of Region IDs. 

###  Visualize polygons before masking
Visually presents the bounding box and polygons to be masked. 

_Any part of a polygon that is outside of the bounding box will be lost. _

### Create Masks
``SET path_to_pkls_folder`` to a folder where you wish to store your regional masks (one .pkl for each polygon region)
   
The big idea:

The given polygons (ex: _Arrigo et al 2017 Sample Areas.shp_) are mapped to the Ocean Colour lat/lon grid and saved (via a pkl file) as an numpy masked array. This runs multiple regions in parallel. 

### Plot
Visually confirm that the stored masks are as expected.



## Step 2: download/retrieve_opendap_gl_chla.ipynb

### Load Regional Masks
In step 1 you created and stored numpy masks for each region as a .pkl file. Now we will load them back into memory. CHANGE path_to_pkls_folder such that it points to the folder containing your pkl files. 


_Note: If your masked files are not sequentially named then update the for x in range(6)_






