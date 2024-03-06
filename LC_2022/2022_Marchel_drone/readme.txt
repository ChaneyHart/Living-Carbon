README document

## desctriptions of python scripts used for multispectral and thrermal data cleaning and extraction

1. cluster_raster.py - used to perform unsupervised classification as an attempt to mask poplar trees automatically (unsuccessful)

2. create_vi_bands.py - for adding 12 vegetation indices (mean and std) to a multispectral raster that has previously been radiometrically corrected (contains reflectance values)

3. extract_temps.py - for extraction of temperature in degrees C (mean and std) within polygons that represent tree boundaries. Takes in a multiband thermal raster where the bands represent separate observations and consolidates the statistics in a csv table.

4. generate_shx_shapefile.py - for checking if a shapefile contains it's associated shx file and attempting to generate a new one if it is missing

5. merge_data.py -  misc. functions for merging and concatenating tabular data using pandas and numpy

6. radiometric_correction.py - takes in raw multispectral and uses empirical line method to convert raw intensity values via images taken of a calibrated reflectance panel before and after the UAS flight mission

7. randomforest_classifier.py - similar to cluster_raster.py, uses randomforest algorithm for unsuccessful unsupervised classification

8. zonal_stats.py - for extraction of mean and standard deviation of pixels within plant boundary polygons for each band of a corrected multispectral with added vegetation bands. Hard-coded for 19 band raster but could easily be modified for more or fewer bands.

9. shadow_area.py  -  for creating bounding boxes and calculating shadowed pixels and shadowed area