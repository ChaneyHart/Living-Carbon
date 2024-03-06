#Extract multispectral from Marchel flight in 2022 for large buffer around each tree (to capture trends in surrounding vegetation)

#Code modified from Cory Garms

#dependencies
library(sf)
library(sp)
library(rgeos)
library(raster)
library(terra)
#install.packages("doParallel")
library(doParallel)
library(stringr)
#install.packages("exactextractr")
library(exactextractr)


#set output filename
outname <- '28AUG2022_Veg_Multispectral_UTM'

#read in points
inpoints = "LC_2022/2022_Marchel_drone/shapefiles/1245_LC_tree_centers_UTM.shp" #shapefile containing tree points
plant_centers  <- st_read(inpoints)

#create circular buffers with 1m radius around points
veg_buffer = st_buffer(plant_centers, dist = 2) #width is in meters

st_layers(dsn = "LC_2022/2022_Marchel_drone/shapefiles/LC_vectors.gpkg")

ROI <- st_read(dsn = "LC_2022/2022_Marchel_drone/shapefiles/LC_vectors.gpkg", layer = "LC_ROI")
ROImask = raster(extent(ROI))
crs(ROImask) = crs(ROI)

# read in raster data
spctrl_file = "LC_2022/2022_Marchel_drone/multispectral/1245_LC_MS_UTM_reflectance.tif"
stackutm <- stack(spctrl_file,bands=c(1:5))#stack(rlist$rlist)
# add in conversion to reflectance?
#add in NDVI and TGI
names(stackutm) = c("b","g","r","ir","re")
sr <- CRS("+proj=utm +zone=10 +ellps=GRS80 +datum=NAD83 +units=m +no_defs") #reference system for projections

stackutm$TGI <- ((stackutm$g)-(0.39*stackutm$r)-(0.61*stackutm$b))#make TGI layer
stackutm$GRVI <- ((stackutm$g - stackutm$r)/(stackutm$g + stackutm$r))#make GRVI layer
stackutm$NDVI <- ((stackutm$ir-stackutm$r)/(stackutm$ir+stackutm$r))
stackutm$NDRE <- ((stackutm$ir-stackutm$re)/(stackutm$ir+stackutm$re))
stackutm$GNDVI <- ((stackutm$ir-stackutm$g)/(stackutm$ir+stackutm$g))

STACK.test = crop(stackutm, extent(ROImask))



print("Beginning zonal statistics calculation")
ptm <- proc.time()
# grab cell number
cell = exact_extract(STACK.test, veg_buffer)
# create a raster with only those cells
no_cores <- detectCores()
cl2 <- makeCluster(no_cores,timeout=180)
registerDoParallel(cl2)
crn_stats = foreach(ii = 1:length(cell),.packages='raster',.combine=rbind,.inorder=T) %do% {
  r = rasterFromCells(STACK.test,ii,values=F)
  foreach(jj = 1:dim(STACK.test)[3],.packages='raster',.combine=cbind,.inorder=T) %dopar% {
    mean(getValues(crop(STACK.test[[jj]],r)),na.rm=TRUE)
  }
}

crn_step1 = foreach(ii = 1:length(cell),.packages='raster',.combine=rbind,.inorder=T) %do% {
  r = rasterFromCells(STACK.test,ii,values=F)
}


crn_step1 = foreach(i = 1:length(cell),.combine=rbind) %do% {
  r = rasterFromCells(STACK.test,i, values = F)
}


r_test <- rasterFromCells(STACK.test,4566,values = F)
mean(getValues(crop(STACK.test,r_test)),na.rm=TRUE)


proc.time() - ptm
stopCluster(cl2)
crn_stats = data.frame(crn_stats)
names(crn_stats) = c("b","g","r","ir","re", "TGI", "GRVI", "NDVI", "NDRE", "GNDVI")
crn_stats$id = c(1:dim(crn_stats))
