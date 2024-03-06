from osgeo import gdal
import numpy as np

input_im = gdal.Open('../Multispectral/LC/1245_LC_MS_UTM.tif')

output_file = '../Multispectral/1245_LC_MS_UTM_reflectance.tif'
driver_tiff = gdal.GetDriverByName('GTiff')

print(dir(input_im))

output_im = driver_tiff.Create(output_file,xsize= input_im.RasterXSize, ysize= input_im.RasterYSize, bands=6, eType= gdal.GDT_Float32)
output_im.SetProjection(input_im.GetProjection())
output_im.SetGeoTransform(input_im.GetGeoTransform())


print(dir(output_im))

print(output_im.GetLayer())

bands = input_im.RasterCount

panelCalibration = { 
    "Blue": 0.62, 
    "Green": 0.46, 
    "Red": 0.61, 
    "Red edge": 0.55, 
    "NIR": 0.36,
    "Alpha" : 1
}

panelEmpirical = {
    "Blue": 40395, 
    "Green": 28834, 
    "Red": 38839, 
    "Red edge": 33324, 
    "NIR": 22028,
    "Alpha" : 1
}

mult_list = []

for key in panelCalibration.keys():
    ref = panelCalibration.get(key)
    emp = panelEmpirical.get(key)
    mult_list.append(ref/emp)

# print(mult_list)


print("[ RASTER BAND COUNT ]: "+ str(bands))
for band in range( bands ):
    band += 1



    # print "[ GETTING BAND ]: ", band
    srcband = input_im.GetRasterBand(band)

 

    stats = srcband.GetStatistics( 0, 1 )

    print("[ STATS ] =  Minimum=%.3f, Maximum=%.3f, Mean=%.3f, StdDev=%.3f" % ( \
            stats[0], stats[1], stats[2], stats[3] ))
    
    input_im_band_ar = srcband.ReadAsArray()
    
    output_im_band_ar = (input_im_band_ar * mult_list[band-1])

    print(output_im_band_ar.max())

    output_im_band_ar = (output_im_band_ar - output_im_band_ar.min()) / (output_im_band_ar.max() - output_im_band_ar.min())

    print(output_im_band_ar.max())

    output_im.GetRasterBand(band).WriteArray(output_im_band_ar)

    newstats = output_im.GetRasterBand(band).GetStatistics(0, 1)

    print("[ STATS ] =  Minimum=%.3f, Maximum=%.3f, Mean=%.3f, StdDev=%.3f" % ( \
        newstats[0], newstats[1], newstats[2], newstats[3] ))

output_im = input_im_band_ar = None
    







