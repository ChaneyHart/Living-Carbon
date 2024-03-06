import geopandas as gpd
import pandas as pd
import numpy as np
import math

#load in tree polygons layer
# shpfile = gpd.read_file('../Shapefiles/PP_mask_5cm_grid.shp')
shp = gpd.read_file('../Shapefiles/LC_mask.shp')
shp_east = gpd.read_file('../Shapefiles/LC_mask_east.shp')

#read rows and columns
rows = gpd.read_file('../Shapefiles/LC_rows.shp')
rows_east = gpd.read_file('../Shapefiles/LC_rows_east.shp').to_crs('EPSG:26910')
cols = gpd.read_file('../Shapefiles/LC_columns.shp')

#read reference csv tabels
r = '../Tables/LC_inventory.csv'
r_east = '../Tables/LC_inventory_east.csv' 

#spatial join using rows and columns
def add_grid(lyr, rows, cols):
    grid_shp = gpd.sjoin_nearest(lyr, rows, how= 'left')
    grid_shp = grid_shp.rename(columns={'fid':'row'})
    grid_shp = grid_shp.drop(columns=['id', 'index_right'])
    grid_shp['row'] = grid_shp['row'].astype(np.int64)
    grid_shp = gpd.sjoin_nearest(grid_shp, cols, how= 'left')
    grid_shp = grid_shp.drop(columns=['index_right', 'fid'])
    grid_shp['column'] = grid_shp['column'].astype(np.int64)

    grid_shp = grid_shp.sort_values(['row', 'column'], ascending=[True, True])

    return grid_shp

#merge gdf with csv data using row and column numbers
def gridmerge(table, reference):
    ref = pd.read_csv(reference)
    merged = pd.merge(table, ref , how='left', on= ['row', 'column'])
    return merged

#merge two csv based on geometry column
def merge_geoms(table, reference):
    tab = pd.read_csv(table)
    ref = pd.read_csv(reference)
    # # these column names are hard coded:
    #for single thermal flight
    ref = ref[['geometry', '1335_degC_mean', '1335_degC_std']]
    #for time series
    # ref = ref[['geometry', '0800_degC_mean', '0800_degC_std', '0900_degC_mean', '0900_degC_std', '1000_degC_mean', '1000_degC_std', '1100_degC_mean', '1100_degC_std', '1200_degC_mean', '1200_degC_std' ]]
    
    merged = pd.merge(tab, ref , how='left', on= ['geometry'])
    print(merged)
    return merged

####THIS SECTION IS FOR MERGING ROWS, COLS, & JOINING WITH OSU DATA#####

# gridlyr = add_grid(shp, rows, cols)
# grid_east = add_grid(shp_east, rows_east, cols)

# merged = gridmerge( gridlyr , r )
# merged_east = gridmerge( grid_east , r_east )

# print(merged)
# print(merged_east)

# merged.to_file('../Shapefiles/LC_bigblock_trees.shp')
# merged.to_file('../Shapefiles/LC_vectors.gpkg', layer='LC_trees_bigblock', driver='GPKG')
# merged.to_csv('../Tables/LC_trees_bigblock.csv')

# merged_east.to_file('../Shapefiles/LC_smallblock_trees.shp')
# merged_east.to_file('../Shapefiles/LC_vectors.gpkg', layer='LC_trees_smallblock', driver='GPKG')
# merged_east.to_csv('../Tables/LC_trees_smallblock.csv')

#########################

#########THIS SECTION IS FOR JOINING THERMAL TIMESERIES WITH MULTISPECTRAL DATA##############
#purple poplar
ms_PP = '../Tables/UM-PP_zonal_stats.csv'
tir_PP = '../Tables/UM-PP_degC.csv'

# # big block
# ms_big = '../Tables/LC_zonal_stats_bigblock.csv'
# tir_big = '../Tables/LC_degC_bigblock.csv'
# # small block
# ms_small = '../Tables/LC_zonal_stats_smallblock.csv'
# tir_small = '../Tables/LC_degC_smallblock.csv'

#merge thermal and ms for big and small 
# ms_tir_big = merge_geoms(ms_big, tir_big)
# ms_tir_small = merge_geoms(ms_small, tir_small)

#merge for pp
ms_tir_PP = merge_geoms(ms_PP, tir_PP)

#fix column name
# ms_tir_small = ms_tir_small.rename(columns={'Event ':'Event'})

# #add block column
# ms_tir_big['block'] = 'big'
# ms_tir_small['block'] = 'small'
ms_tir_PP['block'] = 'PP'

#concatenate big and small blocks
# ms_tir_big = pd.concat([ms_tir_big, ms_tir_small])

#write out to csv
# ms_tir_small.to_csv('../Tables/LC_MS_TIR_smallblock.csv')
# ms_tir_big.to_csv('../Tables/LC_MS_TIR_smallblock.csv')
ms_tir_PP.to_csv('../Tables/UM-PP_MS_TIR_complete.csv')



# ms_tir_merged.to_csv('../Tables/LC_MS_TIR_complete.csv')





