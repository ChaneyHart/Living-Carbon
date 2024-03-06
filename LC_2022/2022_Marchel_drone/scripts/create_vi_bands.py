import numpy as np
import rasterio





def add_veg_indices(input_file, output_file):
    # Open the input TIFF file
    with rasterio.open(input_file) as src:
        # Read the raster data into a numpy array
        data = src.read()

        # Extract the individual bands
        # Extract the red and near-infrared bands
        blue = data[0, :, :]
        green = data[1, :, :]
        red = data[2, :, :]
        ir = data[3, :, :]
        re = data[4, :, :]

        #normalized rgb
        b_norm = blue / (red + green + blue)
        g_norm = green / (red + green + blue)
        r_norm = red / (red + green + blue)

        #calculate and clean ratio arrays
        green_ratio = np.divide(np.ones_like(green), green, out= np.zeros_like(green), where=green!=0)
        re_ratio = np.divide(np.ones_like(green), re, out= np.zeros_like(re), where=re!=0)

        # Calculate the indices
        #RGB indices
        tgi = (green - (0.39*red) - (0.61 * blue))
        grvi = (green- red) / (green + red)
        mgrvi = ((green)**2 - (red)**2) / ((green)**2 + (red)**2)
        exg = 2*g_norm - r_norm - b_norm
        exgr = exg - 1.4*r_norm - g_norm

        #near ir and rededge indices
        ndvi = (ir - red) / (ir + red)
        ndre = (ir - re) / (ir + re)
        gndvi = (ir - green) / (ir + green)
        savi = ((ir - red)*(1 + 0.5))/(ir + red + 0.5)
        osavi = ((ir - red)*(1 + 0.16))/(ir + red + 0.16)
        msavi = ((2*ir + 1) - ((2*ir + 1)**2 - 8*(ir-red))**0.5)/2
        gci = (ir/green) - 1
        reci = (ir/re) - 1
        ari = green_ratio - re_ratio
        
        # Create a new array with 19 bands
        new_image = np.concatenate((data[0:5], np.expand_dims(tgi, axis=0), np.expand_dims(grvi, axis=0), np.expand_dims(mgrvi, axis=0), np.expand_dims(exg, axis=0), np.expand_dims(exgr, axis=0), 
                                      np.expand_dims(ndvi, axis=0), np.expand_dims(ndre, axis=0), np.expand_dims(gndvi, axis=0), np.expand_dims(savi, axis=0), np.expand_dims(osavi, axis=0),
                                        np.expand_dims(msavi, axis=0), np.expand_dims(gci, axis=0), np.expand_dims(reci, axis=0), np.expand_dims(ari, axis=0) ), axis=0)

        new_profile = src.profile
        new_profile.update(count=new_image.shape[0])   
        
        # Create the output TIFF file with 6 bands
        with rasterio.open(output_file, 'w', **new_profile) as dst:
            dst.write(new_image)

# Usage example
# input_file = '../Multispectral/UM/1320_UM_MS_clip.tif'
input_file = '../Multispectral/UM/1320_UM_MS_UTM_reflectance.tif'

output_file = '../Multispectral/UM/1320_UM_MS_full_19band.tif'
add_veg_indices(input_file, output_file)


