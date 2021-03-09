
###############################################################
##  This program calculates LG for SAV and Water line        ##
##  Then it calculates PSI for all data and output to file   ##
##  Then it tests PSI (both) for accuracy in mapping         ##
##  Fields in csv file: first is the unique ID called FID    ##
##  then all the bands with wavelength as the band name      ##
##  Column with classes "WAT" and "SAV" is called "Class"    ##
###############################################################


require("stats")
require("ggplot2")

###############   MIGHT REQUIRE CHANGING    #################

# set working directory
setwd("X:/delta_sav/code/R/ImageAnalysis2020/")

# directory collecting all the files
outdir = "X:/delta_sav/raster/analysis/indx_waterline/202007/"
dirlen = nchar(outdir)

# path and name of csv file with all water and SAV pixel spectra  
csv_filename = paste(outdir, "R_SAVwaterpts2020_bal.csv", sep='')
#csv_trnfile  = paste(outdir, "watersav_201810_all.csv", sep='')
csv_trnfile = csv_filename

# path and filename of all output files
sumFilePSIW = paste(substr(csv_filename, 1, dirlen+17), "_test_PSIW_accuracy.csv", sep='')
sumFilePSIS = paste(substr(csv_filename, 1, dirlen+17), "_test_PSIS_accuracy.csv", sep='')
indxValPSIW = paste(substr(csv_filename, 1, dirlen+17), "_test_PSIW_calcindx.csv", sep='')
indxValPSIS = paste(substr(csv_filename, 1, dirlen+17), "_test_PSIS_calcindx.csv", sep='')
rsqFilePSIW = paste(substr(csv_filename, 1, dirlen+17), "_test_PSIW_rsquared.csv", sep='')
rsqFilePSIS = paste(substr(csv_filename, 1, dirlen+17), "_test_PSIS_rsquared.csv", sep='')

# what are water and SAV classes called?
water = "Wat"
sav   = "SAV"

# ignore the edge pixels of the sorted distribution

ignorePix = 3000

#####################     END CHANGE      ######################

# read file
df.in    <- read.csv(csv_filename, header=T)
records  <- nrow(df.in)
conames  <- colnames(df.in)

# read training file to get regression coefficients
df.train <- read.csv(csv_trnfile, header=T)
rectrn   <- nrow(df.train)
conames  <- colnames(df.train)

# convert to dataframe
# df.in  <- data.frame(infile)
df.wat <- data.frame(subset(df.in, Class == water))
df.sav <- data.frame(subset(df.in, Class == sav))

df.trnwat <- data.frame(subset(df.train, Class == water))
df.trnsav <- data.frame(subset(df.train, Class == sav))

# number of records of water or SAV
recwat <- nrow(df.wat)
recsav <- nrow(df.sav)

# assuming that all SAV records come before water records, stno to enno is 
# the middle of distribution ignoring ignorePix number of pixels at edge
stno = 1 + ignorePix
enno = records - ignorePix
# increment to test for accuracy
incrm = 10

# get column names
conames <- colnames(df.in)
cnmstrn <- colnames(df.train)

# get number of columns
no_cols <- ncol(df.in)
# get column number for "Class" column
classnum = match("Class", conames)

#########################################
# get start and end column numbers for  #
# PSIW: bands 50-57 vs. bands 65-68     #
# PSIS: bands 31-46 vs. bands 65-68     #
# This only applies to AVIRIS-ng files  #
#########################################

# sav1 = match("band31", conames)
# sav2 = match("band46", conames)
# sav3 = match("band65", conames)
# sav4 = match("band68", conames)
# wat1 = match("band50", conames)
# wat2 = match("band57", conames)
# wat4 = match("band68", conames)
# wat3 = sav3

# new assignment for this year - HyMap data
########################################
# get start and end column numbers for #
# PSIW: bands 567-663 vs. 670-753      #
# PSIS: bands 567-663 vs. 670-753      #
# This only applies to HyMap files     #
########################################

sav1 = 5
sav2 = 19
sav3 = 20
sav4 = 32
wat3 = sav3
wat4 = sav4
wat1 = sav1
wat2 = sav2

# column number + stband gives the actual band number
#stband = 31 - sav1
stband = 0

# create an empty matrix to fill up with r squares, p values, and max accuracies of working combinations
Summary_SVIW <- matrix(, nrow = 200, ncol = 7)
Summary_SVIS <- matrix(, nrow = 200, ncol = 7)

# initialize counter for populating above matrices with working combinations
k = 1   # counter for water
l = 1   # counter for sav

# column names for above matrices
bnames = c("band_y", "band_x", "r_sqared", "p_value", "slope", "intercept", "class_acc")
colnames(Summary_SVIW) = bnames
colnames(Summary_SVIS) = bnames

# data frame which will keep increasing as new combinations are tested for SVIW
indxval_wat = as.data.frame(cbind(df.in[,1], as.character(df.in[,classnum])))
indxval_sav = as.data.frame(cbind(df.in[,1], as.character(df.in[,classnum])))
indxcolwat = c("FID", "Class")
indxcolsav = c("FID", "Class")

# create an empty matrix to fill up with r squares for each and every combination
rsq_savline <- matrix(, nrow = (sav2 - sav1 + 1), ncol = (sav4 - sav3 + 1))
rsq_watline <- matrix(, nrow = (wat2 - wat1 + 1), ncol = (wat4 - wat3 + 1))

# column and row names reflect actual band numbers
col_names_sav = conames[sav3:sav4]
row_names_sav = conames[sav1:sav2]
col_names_wat = conames[wat3:wat4]
row_names_wat = conames[wat1:wat2]

# assign column and row names to the rsq matrices
colnames(rsq_savline) = conames[sav3:sav4]
rownames(rsq_savline) = conames[sav1:sav2]
colnames(rsq_watline) = conames[wat3:wat4]
rownames(rsq_watline) = conames[wat1:wat2]

# temporary array to get max accuracy
tempacc  = vector(mode="double", records/incrm)
indxSVIW = vector(mode="double", records)
indxSVIS = vector(mode="double", records)

# Function to calculate the SAV Vertical Index (SVI)
SVI_calculate <- function(x, y, slope, inter) {
  
  SVI <- abs(slope*x - y + inter)/(sqrt(slope*slope + 1))
  return(SVI)
}

# FOR SVIS, do these bands
for (i in sav1:sav2) {
  for (j in sav3:sav4) {

    # generate formula to use in lm function
    formula <- as.formula(paste(conames[i], ' ~ ', conames[j], sep=""))

    col1 = match(conames[i], cnmstrn)
    col2 = match(conames[j], cnmstrn)
    
    # linear regressions between current "i" and "j"
    lrmod_sav <- lm(formula, data = df.trnsav)
    
    coeff.savline = coef(lrmod_sav)
    
    coeff.int = coeff.savline[[1]][1]
    coeff.slp = coeff.savline[[2]][1]
    
    basic.aes = aes(x = df.in[[i]], y = df.in[[j]], colour = Class)
    
    ggplot(df.in, basic.aes) + geom_point(size = 1.5) + theme_bw(base_size = 20) +
      theme(panel.grid.major = element_blank(), legend.position = "none") + 
      scale_color_manual(values = c("#CC0000", "#0066FF")) + geom_smooth(method=lm, se=FALSE) +
      geom_abline(intercept = coeff.int, slope = coeff.slp)
    
    basic.aes = aes(x = df.train[[col1]], y = df.train[[col2]], colour = Class)

    ggplot(df.train, basic.aes) + geom_point(size = 1.5) + theme_bw(base_size = 20) +
      theme(panel.grid.major = element_blank(), legend.position = "none") + 
      scale_color_manual(values = c("#CC0000", "#0066FF")) + geom_smooth(method=lm, se=FALSE)
      #geom_abline(intercept = coeff.int, slope = coeff.slp)
    
    
    # get r-squared value
    rsqsav = summary(lrmod_sav)$r.squared
    rsq_savline[i-sav1+1, j-sav3+1] = rsqsav

    this_colname = paste("SVIS", "_", substr(conames[i], 5, 6), "_", substr(conames[j], 5, 6), sep="")
    
    # if r.squared is > 0.8, then run SVIS (SAV vertical index using sav line)
    if (rsqsav > 0.7) {
      
      # get f-statistic to find out p-value
      fstat <- summary(lrmod_sav)$fstatistic
      
      # get p-value from F distribution
      pval  <- pf(fstat[1],fstat[2],fstat[3],lower.tail=F)
      
      # equation is band_y = coeff[1,1](intercept) + coeff[2,1](slope)*band_x
      coeff <- summary(lrmod_sav)$coefficients
      slope <- coeff[2,1]
      inter <- coeff[1,1]
      
      # note all necessary parameters in comprehensive results table
      Summary_SVIS[l, 1] = i + stband
      Summary_SVIS[l, 2] = j + stband
      Summary_SVIS[l, 3] = rsqsav
      Summary_SVIS[l, 4] = pval
      Summary_SVIS[l, 5] = slope
      Summary_SVIS[l, 6] = inter
      
      # calculate SVIS for each SAV and water pixel in df.in
      indxSVIS <- SVI_calculate(df.in[[j]], df.in[[i]], slope, inter)
      # combine calculated SVIW with class info
      indxval_sav <- cbind(indxval_sav, indxSVIS)
      # assign column names
      indxcolsav <- c(indxcolsav, this_colname)
      colnames(indxval_sav) = indxcolsav
      
      # test accuracy using this SVIS
      
      # SORT the vector by index
      indxSVIS_sort <- indxval_sav[order(indxval_sav[[l+2]]),]

      # intialize accuracy array
      tempacc[] = 0
      
      for (r in seq(stno, enno, incrm)) {
        
        # move threshhold from smallest value to biggest value
        threshold = indxSVIS_sort[r, l+2]
        # to count the number of records of a certain class
        cor_sav = nrow(subset(indxSVIS_sort, (indxSVIS_sort[[2]] == sav)   & (indxSVIS_sort[[l+2]] <= threshold)))
        cor_wat = nrow(subset(indxSVIS_sort, (indxSVIS_sort[[2]] == water) & (indxSVIS_sort[[l+2]] >  threshold)))
        # new accuracy
        tempacc[(r-stno)/incrm + 1] = ((cor_sav + cor_wat)/records)*100

      }
      
      # store maximum accuracy
      Summary_SVIS[l, 7] = max(tempacc)
      
      # advance counter for working combinations populating matrix Summary_SVIS
      l = l + 1
    }   
    cat("completed processing: band1 - ", i+stband, "band2 - ", j+stband, "\n")
  }
}

# Write all relevant information of the PSIS data
write.csv(Summary_SVIS, file=sumFilePSIS)
write.csv(indxval_sav,  file=indxValPSIS)
write.csv(rsq_savline,  file=rsqFilePSIS)

# FOR SVIW, do these bands
for (i in wat1:wat2) {
  for (j in wat3:wat4) {

    # generate formula to use in lm function
    formula <- as.formula(paste(conames[i], ' ~ ', conames[j], sep=""))
    
    # linear regressions between all i's and j'
    lrmod_wat <- lm(formula, data = df.trnwat)
    
    # get r-squared value
    rsqwat = summary(lrmod_wat)$r.squared
    rsq_watline[i-wat1+1, j-wat3+1] = rsqwat
    
    this_colname = paste("SVIW", "_", substr(conames[i], 5, 6), "_", substr(conames[j], 5, 6), sep="")
    
    lrmod_wat
    
    # if r.squared is > 0.8, then run SVIW (SAV vertical index using Water line)
    if (rsqwat > 0.75) {
      
      # get f-statistic to find out p-value
      fstat <- summary(lrmod_wat)$fstatistic
      
      # get p-value from F distribution
      pval  <- pf(fstat[1],fstat[2],fstat[3],lower.tail=F)
      
      # equation is band_y = coeff[1,1](intercept) + coeff[2,1](slope)*band_x
      coeff <- summary(lrmod_wat)$coefficients
      slope <- coeff[2,1]
      inter <- coeff[1,1]
      
      # note all necessary parameters in comprehensive results table
      Summary_SVIW[k, 1] = i + stband
      Summary_SVIW[k, 2] = j + stband
      Summary_SVIW[k, 3] = rsqwat
      Summary_SVIW[k, 4] = pval
      Summary_SVIW[k, 5] = slope
      Summary_SVIW[k, 6] = inter
      
      # calculate SVIW for each SAV and water pixel in df.in
      indxSVIW <- SVI_calculate(df.in[[j]], df.in[[i]], slope, inter)
      # combine calculated SVIW with class info
      indxval_wat <- cbind(indxval_wat, indxSVIW)
      # assign column names
      indxcolwat <- c(indxcolwat, this_colname)
      colnames(indxval_wat) = indxcolwat
      
      # test accuracy using this SVIW
      
      # SORT the vector by index
      indxSVIW_sort <- indxval_wat[order(indxval_wat[[k+2]]),]
      
      # intialize accuracy array
      tempacc[] = 0
      
      for (r in seq(stno, enno, incrm)) {
        
        threshold = indxSVIW_sort[r, k+2]
        # to count the number of records of a certain class
        cor_wat = nrow(subset(indxSVIW_sort, (indxSVIW_sort[[2]] == water) & (indxSVIW_sort[[k+2]] <= threshold)))
        cor_sav = nrow(subset(indxSVIW_sort, (indxSVIW_sort[[2]] == sav)   & (indxSVIW_sort[[k+2]] >  threshold)))
        # new overall accuracy for two classes
        tempacc[(r-stno)/incrm + 1] = ((cor_sav + cor_wat)/records)*100
      }
      
      # store maximum accuracy
      Summary_SVIW[k, 7] = max(tempacc)
      
      # advance counter for working combinations populating matrix Summary_SVIW
      k = k + 1
      
    } # end if conditions for rsq > 0.8
    
    cat("completed processing: band1 - ", i+stband, "band2 - ", j+stband, "\n")
    
  } # end j loop for columns for SVIW
} # end i loop for rows for SVIW


write.csv(Summary_SVIW, file=sumFilePSIW)
write.csv(indxval_wat,  file=indxValPSIW)
write.csv(rsq_watline,  file=rsqFilePSIW)


###############################################################
####  MAKING THE SCATTER-PLOT DEMO FOR SAV AND WATER LINE  ####
###############################################################

# calculate coefficients for the SAV line 

require("ggplot2")

lm.sav = lm(X0.656 ~ X0.698, data = df.sav)

coeff.savline = coef(lm(X0.656 ~ X0.698, data = df.sav))

coeff.int = coeff.savline[[1]][1]
coeff.slp = coeff.savline[[2]][1]

basic.aes = aes(x = X0.656, y = X0.698, colour = Class)

ggplot(df.in, basic.aes) + geom_point(size = 1.5) + theme_bw(base_size = 20) +
      theme(panel.grid.major = element_blank(), legend.position = "none") + 
      scale_color_manual(values = c("#CC0000", "#0066FF")) + geom_smooth(method=lm, se=FALSE)
      
  
#  geom_abline(intercept = coeff.int, slope = coeff.slp)
#  geom_smooth(method=lm, se=FALSE)

ggsave(filename = "ScatterSAVline.png", plot = last_plot(), width = 6, height = 5, dpi = 400)




# calculate co-efficients for the water line 

lm.sav = lm(X0.649 ~ X0.725, data = df.wat)

coeff.savline = coef(lm(X0.649 ~ X0.725, data = df.wat))

coeff.int = coeff.savline[[1]][1]
coeff.slp = coeff.savline[[2]][1]

basic.aes = aes(x = X0.649, y = X0.725, colour = Class)

ggplot(df.in, basic.aes) + geom_point(size = 1.5) + theme_bw(base_size = 20) +
  theme(panel.grid.major = element_blank(), legend.position = "none") + 
  scale_color_manual(values = c("#CC0000", "#0066FF")) + geom_smooth(method=lm, se=FALSE)


#  geom_abline(intercept = coeff.int, slope = coeff.slp)
#  geom_smooth(method=lm, se=FALSE)

ggsave(filename = "ScatterWaterline.png", plot = last_plot(), width = 6, height = 5, dpi = 400)

