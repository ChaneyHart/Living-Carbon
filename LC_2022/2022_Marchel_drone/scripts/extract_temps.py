import geopandas as gpd
import pandas as pd
import rasterstats
import rasterio

# Path to the shapefile
shapefile_path = "../Shapefiles/LC_vectors.gpkg"

# Path to the raster files 
# rasterfile_path = '../Thermal/1335_UM_TIR_PP_clip_rect.tif'
rasterfile_path = '../Thermal/28AUG2022_LC_TIR_TimeSeries_UTM.tif'

# temp_list = ['1335_degC', 'alpha']

temp_list = ['0800_degC', '0900_degC', '1000_degC', '1100_degC', '1200_degC']

def get_temps(shpdir, rasdir, band_names):

    multiband = rasterio.open(rasterfile_path)

    affine = multiband.transform

    print(multiband.count)

    # Read the shapefile using geopandas
    gdf = gpd.read_file(shapefile_path, layer='LC_trees_smallblock')

    for i in range(multiband.count):
        name = band_names[i]
        print(name)

        # Calculate zonal statistics
        stats = rasterstats.zonal_stats(gdf, multiband.read(i+1), stats=["mean", "std"], affine= affine)

        # Create a pandas dataframe from the statistics
        data = pd.DataFrame(stats)

        #rename columns based on band
        data = data.rename(columns={'mean': f'{name}_mean', 'std': f'{name}_std'})

        # Merge the statistics dataframe with the original shapefile
        gdf = pd.concat([gdf, data], axis=1)

    try:
        gdf = gdf.drop(columns=['alpha_mean', 'alpha_std'])
    except:
        print('there is no Alpha band')

    # Print the merged dataframe
    print(gdf)
    gdf.to_file('../Shapefiles/LC_degC_smallblock.shp')
    gdf.to_file('../Shapefiles/LC_vectors.gpkg', layer='LC_TempC_smallblock', driver='GPKG')
    gdf.to_csv('../Tables/LC_degC_smallblock.csv')

get_temps(shapefile_path, rasterfile_path, temp_list)









