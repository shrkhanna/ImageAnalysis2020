#########################################
##                                     ##
##   EXTRACT VECTOR DATA FROM IMAGES   ##
##  Only valid for POLYGON shapefiles  ##
##       and MULTIBAND images          ##
##     example with SpecTIR data       ##
## Can also divide into test/training  ##
##     developed by Shruti Khanna      ##
##      last updated: Apr, 2021        ##
##                                     ##
## make sure that the shapefile has a  ##
## field called ORIG_FID that contains ##
## a unique FID for each vector point  ##
##                                     ##
#########################################

# IMPORTANT information about your vectors
# create a new field called ORIG_FID in your shapefile that is the same as your FID
# the FID number will change through this script but the ORIG_FID should stay the same

# working directory with the R program file
dir_progr = "X:/delta_sav/code/R/ImageAnalysis2020/"

# set working directory
setwd(dir_progr)

# install the necessary packages and load them
# install.packages(c("raster","sp","rgdal","asbio"))
require("raster")
require("sp")
require("rgdal")
require("rgeos")
library("tidyverse")

# name all necessary folders and files

#########  CHANGE  #########

# image directory from where data will be extracted
dir_image = "X:/delta_sav/raster/classification/allinputs/202007/"
# training/test data directory without the last slash front slash
dir_shape = "X:/delta_sav/vector/field_data/2020/Data3"
# output directory for the csv file
dir_out   = "X:/delta_sav/raster/classification/training_test/202007/"
# suffix of images to be processed
imgsuf  = "_all.img"
# name of shapefile without the .shp extension
name_shape  = "ALL2020polygons_"
# if there is a mask_value to be ignored in the image files
mskval = 0
# begin processing at this file
stfile = 1
# end processing at this file
enfile = 0
# divide into test and training percentages
divide_trntst = TRUE
train_frac = 0.6


######## END CHANGE ########

# name of training csv files
if (divide_trntst == FALSE) {
  name_csv = paste(dir_out, "R_", name_shape, "tsttrn.csv", sep="")
} else {
  name_csv.trn = paste(dir_out, "R_", name_shape, "trn.csv", sep="")
  name_csv.tst = paste(dir_out, "R_", name_shape, "tst.csv", sep="")
}

# get list of files; check if file names make sense
img_list <- list.files(dir_image, pattern = imgsuf, full.names=TRUE)
img_list <- subset(img_list, !grepl(".enp", img_list))
no_files <- length(img_list)

if (enfile < stfile) { enfile = no_files }

#enfile = 33

# read metadata
img_attrib = unique(attr(GDALinfo(img_list[1]), 'mdata'))
# separate band number from band name
img_bnames = sapply(strsplit(img_attrib, '='), function(x)x[2])
# get index of band names (because 51 comes after 5)
index = gsub("Band_", "", sapply(strsplit(img_attrib, '='), function(x)x[1]))
# reorder band names by their index to them back in right order
bnames = img_bnames[order(as.integer(index))]
# get number of bands in image
bnumbr = length(bnames)

# read the data shapefile: readOGR(dsn = destination folder, layer = filename)
vector_shape <- readOGR(dsn = dir_shape, layer = name_shape)

# get names of field data columns
fields <- colnames(vector_shape@data)

first = 0

# for (i in seq(enfile, stfile, by=-1)) {
#     print(i)
#   }

# for giving precedence to the flightline to the south
# for (i in seq(enfile, stfile, by=-1)) {

# for giving precedence to the flightline to the north
for (i in seq(stfile, enfile, by=1)) {
    
  # read the input image (give entire path)
  input_image = brick(img_list[i])
  # intialize number of polygons to 0
  no_poly_trn = 0
  
  # if you need to check projections are same or not
  # vector_shape@proj4string
  # input_image@crs
  
  # get subset of points that intersect with the image i
  vint = raster::intersect(vector_shape, as(extent(input_image), 'SpatialPolygons'))
  
  # test if any vector data intersected with the image i
  if (!is.null(vint)) {
    if (first == 0) {
      # first file found where vector data intersected with image
      # turn flag "ON" so that master array can be initialized
      first = i
    }
    # number of polygons intersecting with image i
    no_poly_trn = nrow(vint)
  }
  
  # if more than 1 polygon intersected with image i, then ...
  if (!is.null(vint)) {
    
    # Extract training data as a data frame
    extrn = raster::extract(input_image, vint, cellnumbers=TRUE)
    # Make a vector of no. of pixels extracted per polygon
    freq = t(as.data.frame(lapply(extrn, length)))
    # divide by number of columns to get number of rows
    freq = freq/(bnumbr+1)
    # Total number of pixels extracted from the image
    total = sum(freq)
    
    # convert large list into a single list
    subtrn = NULL
    for (j in 1:no_poly_trn) {
      subtrn = rbind(subtrn, as.data.frame(extrn[[j]]))
    }    

    # add columns to vector data containing number of pixles in each polygon
    vint_freq = as.data.frame(cbind(vint@data, freq))
    # add name of column containing the pixel count
    names(vint_freq) = c(fields, "freq")
    # make this new vector data frame have the same number of rows as the extracted data (duplicate rows for no. of pixels)
    vector = vint_freq[rep(row.names(vint_freq), vint_freq$freq),]
    
    # bind extracted data to the vector data which are now equal length
    subtrnb = cbind(subtrn, vector)
    # give correct column headings to the composite file
    names(subtrnb) = c("FID", bnames, fields, "freq")
    # remove all "0" rows from result data frame
    subtrnb = subset(subtrnb, !((subtrnb[[2]] == mskval) & (subtrnb[[3]] == mskval) & (subtrnb[[4]] == mskval) & (subtrnb[[5]] == mskval)))

    # if this is the first image where vector data intersected, begin master vector
    if (i == first) {
  
      # initialize the master data frames for test and training
      mastrn <- subtrnb
  
      # add vector data to the band names list
      names(mastrn) = c("FID", bnames, fields, "freq")
  
    } else {
      
      # prevent duplication of data where flightlines overlap
      # remove all fids in following file which have been extracted before
      subtrnb = subtrnb[!(subtrnb$ORIG_FID %in% mastrn$ORIG_FID),]
      # append new extracted data to existing data (previously extracted)
      mastrn <- rbind(mastrn, subtrnb)
    }
    
    # status report
    print(i)

  } # end if no polygons condition
}  # end i for loop

if (divide_trntst == FALSE) {
  # write table to csv
  write.csv(mastrn, file=name_csv)
} else {
  # divide test and training data by ORIG_FID and by fraction in train_frac
  # possible commands to try
  # new_df <- df %>% group_by(ID) %>% sample_n(500)
  # train <- mtcars %>% dplyr::sample_frac(.75)
  # test  <- dplyr::anti_join(mtcars, train, by = 'id')
  
  poly.classid = data.frame(cbind(vector_shape@data$Species_1, vector_shape@data$Species_2, vector_shape@data$ORIG_FID))
  colnames(poly.classid) = c('Species1', 'Species2', 'ORIG_FID')

  ########################### CREATE THE "CLASS" COLUMN WITH TARGET SPECIES ##############################
  
  poly.classid$Class = substring(poly.classid$Species1, 1, 3)

  poly.classid$Class = ifelse(poly.classid$Class == 'Rip', 'Riparian', poly.classid$Class)
  poly.classid$Class = ifelse(poly.classid$Class == 'Sha', 'Shadow',   poly.classid$Class)
  poly.classid$Class = ifelse(poly.classid$Class == 'Wat', 'Water',    poly.classid$Class)
  poly.classid$Class = ifelse(poly.classid$Class == 'Soi', 'Soil',     poly.classid$Class)

  poly.classid$Class = ifelse(poly.classid$Species1 == 'Willow-community', 'Riparian', poly.classid$Class)
  poly.classid$Class = ifelse(poly.classid$Species1 == 'EMR-Arundo', 'Arundo',         poly.classid$Class)
  poly.classid$Class = ifelse(poly.classid$Species1 == 'EMR-Phragmites', 'Phragmites', poly.classid$Class)
  poly.classid$Class = ifelse(poly.classid$Species1 == 'FLT-Pennywort', 'Pennywort',   poly.classid$Class)
  poly.classid$Class = ifelse(poly.classid$Species1 == 'FLT-W-hyacinth', 'WHyacinth',  poly.classid$Class)
  poly.classid$Class = ifelse(poly.classid$Species1 == 'FLT-W-primrose', 'WPrimrose',  poly.classid$Class)
  
  poly.classid$Class = ifelse((poly.classid$Species1 == 'EMR-Tule' | poly.classid$Species1 == 'EMR-Cattail' |
                               poly.classid$Species1 == 'EMR'), 'TuleCat', poly.classid$Class)
  
  ########################################################################################################

  poly.classid.sub = poly.classid[(poly.classid$Class == 'FLT' | poly.classid$Class == 'Mud' | poly.classid$Class == 'Pennywort'),]
  poly.classid.fin <- dplyr::anti_join(poly.classid, poly.classid.sub, by = 'ORIG_FID')
  
  paste(unique(poly.classid.fin$Class))
  
  poly.classid.trn <- poly.classid.fin %>% group_by(Class) %>% dplyr::sample_frac(train_frac)
  poly.classid.tst <- dplyr::anti_join(poly.classid.fin, poly.classid.trn, by = 'ORIG_FID')
  poly.classid.tst <- rbind(poly.classid.tst, poly.classid.sub)
    
  poly.classid.trn$ORIG_FID = as.numeric(poly.classid.trn$ORIG_FID)
  poly.classid.tst$ORIG_FID = as.numeric(poly.classid.tst$ORIG_FID)
  
  mastrn.sampled <- left_join(poly.classid.trn, mastrn, by = 'ORIG_FID')
  mastst.sampled <- left_join(poly.classid.tst, mastrn, by = 'ORIG_FID')
  
  # write both tables to csv
  write.csv(mastrn.sampled, file=name_csv.trn)
  write.csv(mastst.sampled, file=name_csv.tst)
}

table(mastrn.sampled$Class)

mastrn.nobal = mastrn.sampled[(mastrn.sampled$Class == 'Arundo' | mastrn.sampled$Class == 'NPV' | mastrn.sampled$Class == 'Phragmites'),]
mastrn.balan = dplyr::anti_join(mastrn.sampled, mastrn.nobal, by = 'Class')
mastrn.balanced <- mastrn.sampled %>% group_by(Class) %>% sample_n(3500)


