
from osgeo import gdal
from osgeo import ogr


import geopandas as gpd

gdal.SetConfigOption('SHAPE_RESTORE_SHX', 'YES')
print(gdal.GetConfigOption('SHAPE_RESTORE_SHX'))

# Open the shapefile
file = '../Shapefiles/PP_mask_5_31.shp'
shapefile = ogr.Open(file)

# shp = gpd.read_file(file)

# print(shp)