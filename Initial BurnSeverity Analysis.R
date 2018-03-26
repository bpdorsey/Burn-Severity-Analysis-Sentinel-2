#############################################
#################################
#PCA burn severity analysis
#much of this code pulled from www.un-spider.org/node/10953/
#2017 - benjamin.dorsey@pc.gc.ca
#################################

#################################
#load needed packages
library(rgdal)			
library(gdalUtils) 
#test gdal install
 gdalinfo()
# gdal_chooseInstallation('JP2OpenJPEG')
 getOption("gdalUtils_gdalPath")[[1]]
 rgdal::gdalDrivers()$name
library(raster)			# Loads the package "raster"
library(rgeos)			# Loads the package "rgeos"
library(rasterVis)		# Loads the package "rasterVis"
library(maptools)		# Loads the package "maptools"
library(RCurl)			# Loads the package "RCurl"
library(devtools)		# Loads the package "devtools"

projUTM <- "+proj=utm +zone=11 +ellps=GRS80 +datum=NAD83 +units=m +no_defs "
p4326 <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs "
epsg3857 <- "+proj=merc +a=6378137 +b=6378137 +lat_ts=0.0 +lon_0=0.0 +x_0=0.0 +y_0=0 +k=1.0 +units=m +nadgrids=@null +wktext  +no_defs";

#################################
#convert a folder of jp2 images to geotiffs 
#################################
setwd("T:/BurnSeverity/data/SentinelImagery/PreFire")
#outdir <- paste(getwd(),"output",sep="/")
df <- list.files(pattern = ".jp2",recursive = T)
#For each 
for(i in 1:length(df)){
  f <-paste(getwd(),df[i],sep="/")
  #if (file.exists(f))  unlink(f) #OR file.remove(f)
  outfile <- paste(getwd(),sub(".jp2",".tif",df[i]),sep="/")
  try(gdalUtils::gdal_translate(df[i],outfile,of="GTiff",output_Raster=TRUE))  #,verbose=TRUE
  }
# End of For each loop
#################################


#################################
#Lets load an analysis extent to clip all grids to 
#################################
#load clip polygon
extent <- readOGR("T:/BurnSeverity/data/extent/mrg_fires.shp") 		# loads the shapefile
projection(extent) <- projUTM #set the current CRS
# if needed change the coordinate system to match the same from the images
#extent_prj<- gdaltransform(extent,t_srs=epsg3857)
#rprojected <- projectRaster(r, crs=p26911))


#Now we have readable files lets make a stack of the bands we need and calc ndvi and NBR then merge the scenes then clip to the anlaysis extent
#################################
getwd()
setwd("T:/BurnSeverity/data/SentinelImagery/PostFire/S2A_MSIL1C_20170927T185131_N0205_R113_T11UMT_20170927T185820.SAFE")
outdir <- paste(getwd(),"output",sep="/")
df <- list.files(pattern = "_B08.tif",recursive = T,include.dirs = T)
df
fnames <- basename(list.files(pattern = "_B08.tif",recursive = T,full.names = FALSE,include.dirs = FALSE))
fnames
#For each 
i = 1
for(i in 1:length(df)){
  b8f <-paste(getwd(),df[i],sep="/")
  b4f <-paste(getwd(),sub("_B08.","_B04.",df[i]),sep="/")
  b12f <-paste(getwd(),sub("_B08.","_B12.",df[i]),sep="/")
  b2f <-paste(getwd(),sub("_B08.","_B02.",df[i]),sep="/")
  b8 <- raster(b8f)
  b4 <- raster(b4f)
  b2 <- raster(b2f)
  #band 12 has a 20 metre spatial resolution, it needs to match the 10m resolution of bands 4 and 8
  b12 <- resample(raster(b12f),b4,method="ngb")
  
  bstack <- stack(b4,b8,b12)
  #calc vegetstackation metrics
  ndvi <- (b8 - b4) / (b8 + b4)
  nbr <- (b8 - b12) / (b8 + b12)
  evi <- (2.5*(b8 - b4) / (b8 + 6*b4 - 7.5*b2 + 1))
  s <- stack(ndvi, nbr, evi)
  
  #Later lets test to see if this is faster
  #fun_ndvi <- function(x) { (x[2]-x[1])/(x[2]+x[1])}
  #nbr <- function(x) { (x[2]-x[1])/(x[2]+x[1])}
  #ndvi <- calc(bstack, fun_ndvi)
  #################
  
  
  ## masks the image to the extent of the shapefile (product to be masked, shapefile used)
  masked <- mask(s,extent) 
  # crops the image to the extent of the shapefile (product to be masked, shapefile used)
  s <- crop(masked,extent) 
  outfile <- paste(outdir,fnames[i],sep="/")
  writeRaster(s, outfile, format="GTiff",overwrite=TRUE)
}




#################################
#merge and clip final rasters
#################################
vinames <- basename(list.files(path =outdir, pattern = ".tif",recursive = T,full.names = FALSE,include.dirs = FALSE))
for (i in 1:length(vinames)){
  if (i == 1){
    lastrast <- raster(vinames[i])
  } 
  lastrast <- merge(lastrast,raster(vinames[i])
}
lastrast<- crop(lastrast,extent)
outfile <- paste(outdir,"VegMetrics",sep="/")
writeRaster(lastrast, outfile, format="raster",overwrite=TRUE)
#################################



#################################
#Subtract pre fire values from post fire values
#################################
# dNBR is calculated through the difference of pre and post-fire NBR
dNBR <- (pre_fire_NBR) - (post_fire_NBR)




#################################
#reclassify each final raster and vectorize
#################################

#reclass all veg as one value
vegc <- reclassify(veg, c(-Inf,0.2,0, 0.2,0.3,1, 0.3,0.4,2, 0.4,0.5,3, 0.5, Inf, 4))
plot(vegc, col = rev(terrain.colors(4)), main = 'NDVI based thresholding')


# sets the ranges that will be used to classify dNBR information about the ranges used, please see the NBR section on the recommended practice
NBR_ranges <- c(-Inf, -500, -1, -500, -251, 1, -251, -101, 2, -101, 99, 3, 99, 269, 4, 269, 439, 5, 439, 659, 6, 659, 1300, 7, 1300, +Inf, -1) 
# sets a classification matrix
class.matrix <- matrix(NBR_ranges, ncol = 3, byrow = TRUE)
# classification matrix is used to classify dNBR_scaled
dNBR_reclass <- reclassify(dNBR_scaled, class.matrix, right=NA)









#################################
#################################

par(mfrow=c(1,3), mar=c(4, 4, 4, 4))
plot(ndvi)
plot(nbr)
plot(evi)





#####################################  
#####################################  
#####################################  







#####################################  
#LAND SAT 8 imagery
#####################################  
#NDVI
# For Sentinel NIR = 7, red = 3.
ndvi <- NDVI(ss, 3, 7)
plot(ndvi, col = rev(terrain.colors(30)), main = 'NDVI from Sentinel')
# calculates pre-fire NBR
pre_fire_NDVI<- NDVI()
#generate raster of classified ndvi image
veg <- calc(ndvi, function(x){x[x < 0.4] <- NA; return(x)})
plot(veg, main = 'Veg cover')
#reclass all veg as one value
vegc <- reclassify(veg, c(-Inf,0.2,0, 0.2,0.3,1, 0.3,0.4,2, 0.4,0.5,3, 0.5, Inf, 4))
plot(vegc, col = rev(terrain.colors(4)), main = 'NDVI based thresholding')





#####################################  
#NBR
 # (reflectance.jan08NIR-reflectance.jan08SWIR2)/reflectance.jan08NIR+reflectance.jan08SWIR2)

# calculates post-fire NBR
#post_fire_NBR<-(reflectance.feb25NIR-reflectance.feb25SWIR2)/(reflectance.feb25NIR+reflectance.feb25SWIR2)

# dNBR is calculated through the difference of pre and post-fire NBR
dNBR <- (pre_fire_NBR) - (post_fire_NBR)
# masks the image to the extent of the shapefile (product to be masked, shapefile used)
dNBR_masked <- mask(dNBR,extent_prj ) 
# crops the image to the extent of the shapefile (product to be masked, shapefile used)
dNBR_cropped <- crop(dNBR_masked,extent_prj ) 
# plots the cropped image
plot(dNBR_cropped) 

###########################################
# scales the dNBR map by 10^3
#dNBR_scaled <- 1000*dNBR_cropped 
# sets the ranges that will be used to classify dNBR information about the ranges used, please see the NBR section on the recommended practice
NBR_ranges <- c(-Inf, -500, -1, -500, -251, 1, -251, -101, 2, -101, 99, 3, 99, 269, 4, 269, 439, 5, 439, 659, 6, 659, 1300, 7, 1300, +Inf, -1) 
# sets a classification matrix
class.matrix <- matrix(NBR_ranges, ncol = 3, byrow = TRUE)
# classification matrix is used to classify dNBR_scaled
dNBR_reclass <- reclassify(dNBR_scaled, class.matrix, right=NA)
#########################################
#############################################

#########################################
#out put and plotting
#############################################


# builds the attribute table for the legend 
dNBR_reclass <- ratify(dNBR_reclass) 
rat <- levels(dNBR_reclass)[[1]]
# creates the text that will be on the legend
rat$legend  <- c("NA", "Enhanced Regrowth, High", "Enhanced Regrowth, Low", "Unburned", "Low Severity", "Moderate-low Severity", "Moderate-high Severity", "High Severity") 
levels(dNBR_reclass) <- rat 

###############################################
# setting the colors for the severity map
my_col=c("white", "darkolivegreen","darkolivegreen4","limegreen", "yellow2", "orange2", "red", "purple")

# plots the burn severity map with title in two lines; first line: 'Burn Severity'; second line: 'Empedrado, Chile'.
plot(dNBR_reclass,col=my_col,legend=F,box=F,axes=F, main="Burn Severity \nEmpedrado, Chile") 

# plots the legend on the right side of the burn severity map
legend(x='right', legend =rat$legend, fill = my_col, y='right') 
#########################################
#############################################