from osgeo import ogr
import numpy as np
from scipy.spatial.distance import cdist
import geopandas as gpd
from shapely.geometry import Point
import math 
import rasterio
import rasterstats
import pandas as pd

#function to add nearest neigbor to points layer
# shp  = '../Shapefiles/LC_bigblock_shadow_centers.shp'
shp  = '../Shapefiles/LC_smallblock_shadow_centers.shp'
# shp  = '../Shapefiles/UM_PP_shadow_centers.shp'



gdf = gpd.read_file(shp)

ds = ogr.Open(shp)
layer = ds.GetLayer()
# get coordinates of points as 2d array
coords = np.array([(feature.GetGeometryRef().GetX(), feature.GetGeometryRef().GetY()) for feature in layer])
# calculate distance matrix between points
distances = cdist(coords, coords)
# sort distances to each point
sorted_distances = np.argsort(distances, axis=1)
# get 2nd column from sorted matrix (it contains 2nd closest points because 1st is always 0!!!)
closest_points = coords[sorted_distances[:,1]]

origins = []

for i in gdf['geometry']:
    origins.append([i.x, i.y])

nn_dists = []

for index, d in enumerate(origins):
    nn_dists.append(math.dist(d, closest_points[index]))

# print(nn_dists)
# print(len(nn_dists))

gdf['nn_dists'] = nn_dists

print(gdf)

new_gdf = gdf 

# for each row compute the buffer based on nearest neighbor distance
def buffer(row):
     return row.geometry.buffer((row.nn_dists/2), cap_style=3)

new_gdf['geometry'] = new_gdf.apply(buffer, axis=1)

# new_gdf.to_file('../Shapefiles/LC_bigblock_shadow_buffers.shp')
new_gdf.to_file('../Shapefiles/LC_smallblock_shadow_buffers.shp')
# new_gdf.to_file('../Shapefiles/UM_PP_shadow_buffers.shp')



# multiband = rasterio.open('../Multispectral/classification/LC_bigblock_shadow_mask.tif')
multiband = rasterio.open('../Multispectral/classification/LC_smallblock_shadow_mask.tif')
# multiband = rasterio.open('../Multispectral/classification/PP_nir_threshold_shadow_mask.tif')


affine = multiband.transform

stats = rasterstats.zonal_stats(new_gdf, multiband.read(1), stats=["sum"], affine= affine)

sums = pd.DataFrame(stats)

new_gdf['shadow_pix'] = sums['sum']

pixel_area = affine[0] * affine[4]

areas = abs(sums * pixel_area)

new_gdf['shadow_area'] = areas

# new_gdf.to_file('../Shapefiles/LC_bigblock_shadow_areas.shp')
new_gdf.to_file('../Shapefiles/LC_smallblock_shadow_areas.shp')
# new_gdf.to_file('../Shapefiles/UM_PP_shadow_areas.shp')

# new_gdf.to_csv('../Shapefiles/UM_PP_shadow_areas.csv')
# new_gdf.to_csv('../Shapefiles/LC_bigblock_shadow_areas.csv')
new_gdf.to_csv('../Shapefiles/LC_smallblock_shadow_areas.csv')





