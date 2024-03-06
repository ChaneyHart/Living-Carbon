import numpy as np
from sklearn import cluster
from osgeo import gdal, gdal_array
import matplotlib.pyplot as plt
from numpy import inf


# Tell GDAL to throw Python exceptions, and register all drivers
gdal.UseExceptions()
gdal.AllRegister()

# Read in raster image
img_ds = gdal.Open('../Multispectral/1320_UM_MS_UTM_VIs_PP_clip.tif', gdal.GA_ReadOnly)

# band = img_ds.GetRasterBand(11)

# img = band.ReadAsArray()

# img[img < -9999] = -9999

# print(img.max())
# print(img.min())


# X = img.reshape((-1,1))

img = np.zeros((img_ds.RasterYSize, img_ds.RasterXSize, img_ds.RasterCount),
               gdal_array.GDALTypeCodeToNumericTypeCode(img_ds.GetRasterBand(1).DataType))

print(img.shape[2])



for b in range(img.shape[2]):
    img[:, :, b] = img_ds.GetRasterBand(b + 1).ReadAsArray()

new_shape = (img.shape[0] * img.shape[1], img.shape[2])

X = img[:, :, :11].reshape(new_shape)

X[X < -9999] = -9999

k_means = cluster.KMeans(n_clusters=5)
k_means.fit(X)

X_cluster = k_means.labels_
#X_cluster = X_cluster.reshape(img.shape)
X_cluster = X_cluster.reshape(img[:, :, 0].shape)
print(X_cluster.shape)

plt.figure(figsize=(20,20))
plt.imshow(X_cluster, cmap="hsv")

plt.show()

ds = gdal.Open("../Multispectral/1320_UM_MS_UTM_VIs_PP_clip.tif")
band = ds.GetRasterBand(1)
arr = band.ReadAsArray()
[cols, rows] = arr.shape

format = "GTiff"
driver = gdal.GetDriverByName(format)


outDataRaster = driver.Create("./k_means_output_5class.tif", rows, cols, 1, gdal.GDT_Byte)
outDataRaster.SetGeoTransform(ds.GetGeoTransform())##sets same geotransform as input
outDataRaster.SetProjection(ds.GetProjection())##sets same projection as input


outDataRaster.GetRasterBand(1).WriteArray(X_cluster)

outDataRaster.FlushCache() ## remove from memory
del outDataRaster ## delete the data (not the actual geotiff)