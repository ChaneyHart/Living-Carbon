import geopandas as gpd
import pandas as pd
import rasterstats
import rasterio

# Path to the shapefile
shapefile_path = "../Shapefiles/UM_vectors.gpkg"


# Path to the raster file (replace 'your_raster.tif' with your actual raster file name)
# rasterfile_path = '../Multispectral/LC/1245_LC_MS_19band.tif'
rasterfile_path = '../Multispectral/UM/1320_UM_MS_clip_19band.tif'



bandnames = ['blue', 'green', 'red', 'nir', 're', 'TGI', 'GRVI', 'MGRVI', 'EXG', 'EXGR',
             'NDVI', 'NDRE', 'GNDVI', 'SAVI', 'OSAVI', 'MSAVI', 'GCI', 'RECI', 'ARI' ]

def multiband_zstats(shpdir, rasdir, band_names):

    multiband = rasterio.open(rasterfile_path)

    affine = multiband.transform

    # Read the shapefile using geopandas
    # gdf = gpd.read_file(shapefile_path, layer='LC_trees_smallblock')
    gdf = gpd.read_file(shapefile_path, layer='PP_trees')


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

    # Print the merged dataframe
    print(gdf)
    # gdf.to_file('../Shapefiles/LC_zonal_stats_smallblock.shp')
    # gdf.to_file('../Shapefiles/LC_vectors.gpkg', layer='LC_MS_stats_smallblock', driver='GPKG')
    # gdf.to_csv('../Tables/LC_zonal_stats_smallblock.csv')

    gdf.to_file('../Shapefiles/UM-PP_zonal_stats.shp')
    gdf.to_file('../Shapefiles/LC_vectors.gpkg', layer='PP_zonal_stats', driver='GPKG')
    gdf.to_csv('../Tables/UM-PP_zonal_stats.csv')

multiband_zstats(shapefile_path, rasterfile_path, bandnames)









