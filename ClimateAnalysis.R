#########################################################
## ClimateAnalyses.R
## Purpose: Analyze climate change over space (2000-4000m) and time (1978-2014)
## Corresponding author:
## William K. Petry (ETH Zürich, william.petry@usys.ethz.ch)
#########################################################
## README ----
# This script is divided into the following sections:
#
# 1. Preliminaries
#   1.1. Load required packages & functions
#   1.2. Set working directory and data options
# 2. Temperature over space
#   2.1. Download data
#   2.2. Read all rasters into workspace, crop to focal region, and calculate tmean across months
#   2.3. Extract spatially paired values from tmean and elevation raster
#   2.4. Run regression analysis to test for change in tmean over elevation
#   2.5. Plot
# 3. Temperature over time
#   3.1. Download, import, clean data
#   3.2. Subset to stations within focal time period and focal region
#   3.3. Convert data to whole units (default is °C x 100)
#   3.4. Calculate mean of climate variable across focal months
#   3.5. Run regression analysis to test for change in tmean over time
#   3.6. Plot mean temperature change over time
# 4. Precipitation over space
#   4.1. Download, import, clean data
#   4.2. Read in all rasters, crop to focal region, and calculate precip across months
#   4.3. Extract spatially paired values from tmean and elevation raster
#   4.4. Run regression analysis to test for change in precip over elevation
#   4.5. Plot
# 5. Precipitation over time
#   5.1. Download, import, clean data
#   5.2. Subset to stations within focal time period and focal region
#   5.3. Convert data to whole units (default is mm * 10)
#   5.4. Calculate mean of climate variable across focal months
#   5.5. Run regression analysis to test for change in precip over time
#   5.6. Plot
# 6. Snow melt over space
#   6.1. Import data
#   6.2. Run regression analysis to test for change in snow melt over elevation
#   6.3. Plot
# 7. Snow melt over time
#   7.1. Import data
#   7.2. Run regression analysis to test for change in snow melt over time
#   7.3. Plot
# 8. Soil moisture over space
#   8.1. Download, import, clean data
#   8.2. Subset to stations within focal time period & average across growing season
#   8.3. Run regression analysis to test for change in soil moisture over elevation
#   8.4. Plot
# 9. Soil moisture over time
#   9.1. Download data
#   9.2. Subset to stations within focal time period and focal region
#   9.3. Convert moisture data to volumetric percentage
#   9.4. Calculate mean soil moisture across growing season
#   9.5. Run regression analysis to test for change in soil moisture over time
#   9.6. Plot
# 10. Calculate paces of climate change
#   10.1. Temperature pace
#   10.2. Precipitation pace
#   10.3. Snow melt pace
#   10.4. Soil moisture pace
# 
# The folowing data must be downloaded manually into the working directory:
# (1) PRISM 4km digital elevation model grid, .BIL format (keep zipped)
# (2) PRISM 800m digital elevation model grid, .BIL format (keep zipped)
# These files were not available by the PRISM FTP or web service. They may be downloaded
# from PRISM at http://www.prism.oregonstate.edu/normals/
#########################################################
## 1. Preliminaries ----
#########################################################
## All subsequent sections depend on this section being run first. Sections 2-9 depend only on this section. Section 10 requires all previous sections to be run.
## 1.1. Load required packages  & functions ====
library(raster)
library(ggplot2)
library(effects) 
library(car)
library(propagate)
library(dplyr)
library(tidyr)
library(magrittr)

# Only one function ('partial.R2') is needed from package 'asbio' v. 1.3-1
# Define partial.R2 here to avoid external dependencies required by 'asbio'
partial.R2 <- function(nested.lm, ref.lm){
    a <- anova(nested.lm)
    b <- anova(ref.lm)
    length.ref <- length(attributes(ref.lm$terms)$"dataClasses")
    length.nested <- length(attributes(nested.lm$terms)$"dataClasses")
    if(length.nested > length.ref) stop("Specify nested model first in arguements")
    if(length.ref - length.nested > 1) stop("Reference and nested model should only differ with repsect to the presence/absence of one predictor")
    SSE.wo <- tail(a$"Sum Sq", 1)
    SSE.with <- tail(b$"Sum Sq", 1)
    P.R2 <- (SSE.wo-SSE.with)/SSE.wo
    P.R2
}

## 1.2. Set working directory and data options ====
# Create folder for climate data, defaulting to system root
dnld.dir <- "~/PetryDryadData/"  # Change this file path to relocate data download folder
dir.create(dnld.dir)
setwd(dnld.dir)

# Extend URL timeout interval to accommodate large files
options(timeout = 180)

# Specify the focal region
center <- c(-106.98886, 38.95775) # center of focal region (the Rocky Mountain Biological Laboratory)
region.size <- c(1, 1) # 1°x1° focal region
focal.region <- extent(c(center[1]-region.size[1]/2, center[1]+region.size[1]/2,
                         center[2]-region.size[2]/2, center[2]+region.size[2]/2))
region.mult <- 2 # multiplier to double focal region size to 2°x2° when weather stations are scarce

# Specify the focal time period
focal.years <- c(1978, 2014) # start year, end year (numeric, YYYY)
focal.months <- 6:8 # start month:end month (numeric without leading zero)

#########################################################
## 2. Temperature over space ----
#########################################################
## Dataset: PRISM 800m Normals 1981-2010
## Variable: "tmean", mean monthly temperature
## URL: http://prism.oregonstate.edu/documents/PRISM_downloads_FTP.pdf

## 2.1 Download data ====
# Specify, retrieve, and unzip temperature rasters from PRISM FTP
# Automatically skips download if files already exist in the correct folder
if(!file.exists(paste0(dnld.dir, "Temperature/PRISM normals 800m tmean"))){
    dir.create(paste0(dnld.dir, "Temperature/PRISM normals 800m tmean"), recursive = T)
}
setwd(paste0(dnld.dir, "Temperature/PRISM normals 800m tmean"))
TmeanSp.files <- vector(mode = "character")
urls <- vector(mode = "character")
for(i in 1:length(focal.months)){
    # Specify file names
    TmeanSp.files[i] <- paste0("PRISM_tmean_30yr_normal_800mM2_",
                               ifelse(focal.months[i] < 10, paste0(0, focal.months[i]), focal.months[i]),
                               "_bil.zip")
    # Download files and unzip only if file doesn't already exist
    if(!file.exists(paste0(substr(TmeanSp.files[i], 1, nchar(TmeanSp.files[i])-3), "bil"))){
        download.file(paste0("ftp://anonymous@prism.nacse.org:21/normals_800m/tmean/",
                             TmeanSp.files[[i]]), TmeanSp.files[i], mode = "wb")
        unzip(TmeanSp.files[i])
        invisible(file.remove(TmeanSp.files[i]))
    }else{
        message(paste("PRISM normals 800m tmean Raster for month", focal.months[i],
                      "has previously been downloaded"))
    }
}

# NOTE: The digital elevation model (DEM) grid used to estimate the elevation of each raster cell must be downloaded manually from PRISM. Visit http://www.prism.oregonstate.edu/normals/, download the elevation grid in .BIL format at 4km and 800m spatial resolution.

# Transfer manually downloaded DEM files to working directory using interactive dialogue
# This section will prompt the user to select the .ZIP file containing the 4km resolution DEM, then the .ZIP file containing the 800m resulution DEM.
if(!file.exists(paste0(dnld.dir, "Elevation/PRISM 4km DEM/PRISM_us_dem_4km_bil.bil"))){
    print("Choose the 4km resolution DEM from the download location")
    dir.create(paste0(dnld.dir, "Elevation/PRISM 4km DEM"), recursive = T)
    unzip(file.choose(), exdir = paste0(dnld.dir, "Elevation/PRISM 4km DEM/"))
}
if(!file.exists(paste(dnld.dir, "Elevation/PRISM 800m DEM/PRISM_us_dem_800m_bil.bil", sep = "/"))){
    print("Choose the 800m resolution DEM from the download location")
    dir.create(paste0(dnld.dir, "Elevation/PRISM 800m DEM"), recursive = T)
    unzip(file.choose(), exdir = paste0(dnld.dir, "Elevation/PRISM 800m DEM/"))
}

## 2.2. Read all rasters into workspace, crop to focal region, and calculate tmean across months ====
TmeanSp.files <- as.list(gsub(".zip", ".bil", TmeanSp.files))
TmeanSp.rasters <- stack(lapply(TmeanSp.files, function(x){crop(raster(x), focal.region)}))
TmeanSp.mean <- overlay(TmeanSp.rasters, fun = mean)

if(!exists("elev800m.raster")){
    elev800m.raster <- crop(raster(paste0(dnld.dir, "Elevation/PRISM 800m DEM/PRISM_us_dem_800m_bil.bil")), focal.region) 
}

## 2.3. Extract spatially paired values from tmean and elevation raster ====
TmeanSp.df <- data.frame(TmeanSp.mean = TmeanSp.mean@data@values, elev = elev800m.raster@data@values)

## 2.4. Run regression analysis to test for change in tmean over elevation ====
TmeanSp.mod <- lm(TmeanSp.mean~elev, data = TmeanSp.df)
summary(TmeanSp.mod)

## 2.5. Plot ====
gg.TmeanSp <- ggplot(TmeanSp.df, aes(elev, TmeanSp.mean))+
    geom_point(size = 8, shape = 16, color = "grey50")+
    scale_x_continuous(name = "Elevation (m)", limits = c(1800, 4250),
                       breaks = c(2000, 3000, 4000), expand = c(0, 0))+
    scale_y_continuous(name = "Mean temp. (°C)", limits = c(0, 21), breaks = c(0, 5, 10, 15, 20), 
                       expand = c(0, 0))+
    geom_smooth(method = 'lm', size = 3, se = F) +
    theme_classic()
gg.TmeanSp

#########################################################
## 3. Temperature over time ----
#########################################################
## Dataset: NOAA NCDC GHCN v3 Monthly Summaries
## Variable: Mean monthly temperature
## URL: http://www.ncdc.noaa.gov/ghcnm/v3.php
## NOTE: SNOTEL dataset has serious problems with temperature readings leading to a large
## overestimate of the pace of climate warming (see http://dx.doi.org/10.1002/2014GL062803).

## 3.1. Download, import, clean data ====
# Create directories and download documentation
if(!file.exists(paste0(dnld.dir, "Temperature/NOAA GHCN tmean"))){
    dir.create(paste0(dnld.dir, "Temperature/NOAA GHCN tmean"), recursive=T)
    download.file("http://www1.ncdc.noaa.gov/pub/data/ghcn/v3/README",
                  paste0(dnld.dir, "Temperature/NOAA GHCN tmean/README.txt"), quiet=T)
    download.file("http://www1.ncdc.noaa.gov/pub/data/cdo/documentation/GHCNDMS_documentation.pdf",
                  paste0(dnld.dir, "Temperature/NOAA GHCN tmean/GHCNDMS_documentation.pdf"), quiet=T) 
}
if(!file.exists(paste0(dnld.dir, "Temperature/USDA NRCS WCC SNOTEL"))){
    dir.create(paste0(dnld.dir, "Temperature/USDA NRCS WCC SNOTEL"), recursive=T)
    download.file("http://www.wcc.nrcs.usda.gov/snotel/SNOTEL_brochure.pdf",
                  paste0(dnld.dir, "Temperature/USDA NRCS WCC SNOTEL/SNOTEL_brochure.pdf"), quiet=T)
}

# Download, import, and clean monthly data + metadata
if(length(list.files(paste0(dnld.dir, "Temperature/NOAA GHCN tmean"),
                     pattern="\\.qca\\.[a-z]{3}$"))==0){
    download.file("ftp://anonymous@ftp.ncdc.noaa.gov/pub/data/ghcn/v3/ghcnm.tavg.latest.qca.tar.gz",
                  paste0(dnld.dir, "Temperature/NOAA GHCN tmean/ghcnm.tavg.latest.qca.tar.gz"),
                  mode="w")
    setwd(paste0(dnld.dir, "Temperature/NOAA GHCN tmean/"))
    ghcn.files <- untar(tarfile=paste0(dnld.dir,
                                       "Temperature/NOAA GHCN tmean/ghcnm.tavg.latest.qca.tar.gz"),
                        list=T)
    untar(tarfile=paste0(dnld.dir, "Temperature/NOAA GHCN tmean/ghcnm.tavg.latest.qca.tar.gz"),
          compressed="gzip")
    file.rename(from=paste0(dnld.dir, "Temperature/NOAA GHCN tmean/",
                            substr(ghcn.files, 3, nchar(ghcn.files))),
                to=paste0(dnld.dir, "Temperature/NOAA GHCN tmean/",
                          substr(ghcn.files, 25, nchar(ghcn.files))))
    unlink(c(paste0(dnld.dir, "Temperature/NOAA GHCN tmean/", substr(ghcn.files, 3, 24)),
             paste0(dnld.dir, "Temperature/NOAA GHCN tmean/ghcnm.tavg.latest.qca.tar.gz")), recursive=T)
}else{
    ghcn.files <- list.files(paste0(dnld.dir, "Temperature/NOAA GHCN tmean/"),
                             pattern="\\.dat$|\\.inv")
}

ghcn.df <- read.fwf(file=paste0(dnld.dir, "Temperature/NOAA GHCN tmean/",
                                substr(ghcn.files[1], nchar(ghcn.files[1])-33, nchar(ghcn.files[1]))),
                    widths=c(11, 4, 4, rep(c(5, 1, 1, 1), times=12)),
                    col.names=c("ID", "YEAR", "ELEMENT",
                                paste0(rep(c("VALUE", "DMFLAG", "QCFLAG", "DSFLAG"), times=12),
                                       rep(1:12, each=4))),
                    header=F, stringsAsFactors=F)
ghcn.meta <- read.fwf(file=paste0(dnld.dir, "Temperature/NOAA GHCN tmean/",
                                  substr(ghcn.files[2], nchar(ghcn.files[2])-33, nchar(ghcn.files[2]))),
                      widths=c(11, 9, 10, 7, 31, 5, 2, 4, 2, 2, 2, 2, 1, 2, 16, 1),
                      col.names=c("ID", "LATITUDE", "LONGITUDE", "STNELEV", "NAME", "GRELEV", "POPCLS", 
                                  "POPSIZ", "TOPO", "STVEG", "STLOC", "OCNDIS", "AIRSTN", "TOWNDIS", 
                                  "GRVEG", "POPCSS"), 
                      header=F, stringsAsFactors=F, strip.white=T, comment.char="~")
# Clean data: replace missing values with NA
ghcn.df[ghcn.df=="-9999"] <- NA

## 3.2. Subset to stations within focal time period and focal region ====
# To focal time period
ghcn.df.crop <- subset(ghcn.df, subset=YEAR%in%focal.years[1]:focal.years[2], 
                       select=c("ID", "YEAR", "ELEMENT", 
                                paste0(rep(c("VALUE", "DMFLAG", "QCFLAG", "DSFLAG"),
                                           times=length(focal.months)), rep(focal.months, each=4))))

# To focal region (exclude stations with coordinates outside of focal region)
ghcn.spdf <- with(ghcn.meta, SpatialPointsDataFrame(coords=data.frame(LONGITUDE, LATITUDE), 
                                                    data=ghcn.meta[, c(1, 4:16)]))
ghcn.df.focal <- ghcn.df.crop[ghcn.df.crop$ID%in%
                                  ghcn.spdf[!is.na(ghcn.spdf%over%as(focal.region*region.mult, 
                                                                     "SpatialPolygons")), ]@data[, 1], ]
ghcn.spdf.focal <- ghcn.spdf[ghcn.spdf$ID%in%
                                 ghcn.spdf[!is.na(ghcn.spdf%over%as(focal.region*region.mult, 
                                                                    "SpatialPolygons")), ]@data[, 1], ]
ghcn.df.focal$ID <- factor(ghcn.df.focal$ID)
# Print metadata for stations within focal region
ghcn.meta[ghcn.meta$ID%in%levels(ghcn.df.focal$ID), 1:5]

## 3.3. Convert data to whole units (default is °C * 100) ====
cols <- paste0(rep("VALUE", times=length(focal.months)), focal.months)
ghcn.df.focal[cols] <- ghcn.df.focal[cols]/100.0

## 3.4. Calculate mean of climate variable across focal months ====
ghcn.df.focal$MEAN.VALUE <- rowMeans(sapply(as.list(paste0("VALUE", focal.months)), 
                                            get, pos=ghcn.df.focal), 
                                     na.rm=F)

## 3.5. Run regression analysis to test for change in tmean over time ====
ghcn.summary.lm <- lm(MEAN.VALUE~YEAR+ID, data=ghcn.df.focal)
summary(ghcn.summary.lm)
Anova(ghcn.summary.lm, type="III")
partial.R2(update(ghcn.summary.lm, ~. -YEAR), ghcn.summary.lm) # estimate partial R2 accounting for station
ghcn.summary.resids <- with(effect("YEAR", ghcn.summary.lm, partial.residuals=T),  # set up plot data
                            data.frame(YEAR=data$YEAR, RESID=residuals, row.names=NULL))

## 3.6. Plot ====
gg.TmeanTime <- ggplot(ghcn.summary.resids, aes(x=YEAR, y=RESID))+
    geom_point(size=8, shape=16, color="grey50")+
    geom_smooth(method="lm", size=3, se=F, color="dodgerblue3")+
    scale_x_continuous(name="Year", limits=focal.years, breaks=c(1980, 1990, 2000, 2010), expand=c(0, 0))+
    scale_y_continuous(name="Mean temp. (°C)", limits=c(0, 21), breaks=c(0, 5, 10, 15, 20), expand=c(0, 0))+
    theme_classic()
gg.TmeanTime

#########################################################
## 4. Precipitation over space ----
#########################################################
## Dataset: PRISM 800m Normals 1981-2010
## Variable: "precip" monthly total precipiation
## URL: http://prism.oregonstate.edu/documents/PRISM_downloads_FTP.pdf

## 4.1. Download, import, clean data ====
if(!file.exists(paste0(dnld.dir, "Precipitation/PRISM normals 800m ppt"))){
    dir.create(paste0(dnld.dir, "Precipitation/PRISM normals 800m ppt"), recursive=T)
}
setwd(paste0(dnld.dir, "Precipitation/PRISM normals 800m ppt"))
PrecipSp.files <- vector(mode="character")
urls <- vector(mode="character")
for(i in 1:length(focal.months)){
    # Specify file names
    PrecipSp.files[i] <- paste0("PRISM_ppt_30yr_normal_800mM2_", 
                                ifelse(focal.months[i]<10, paste0(0, focal.months[i]), focal.months[i]), 
                                "_bil.zip")
    # Download files and unzip only if file doesn't already exist
    if(!file.exists(paste0(substr(PrecipSp.files[i], 1, nchar(PrecipSp.files[i])-3), "bil"))){
        download.file(paste0("ftp://anonymous@prism.nacse.org:21/normals_800m/ppt/", 
                             PrecipSp.files[[i]]), PrecipSp.files[i], mode="wb")
        unzip(PrecipSp.files[i])
        invisible(file.remove(PrecipSp.files[i]))
    }else{
        print(paste("PRISM normals 800m precip Raster for month", focal.months[i], 
                    "has previously been downloaded"))
    }
}

# NOTE: The digital elevation model (DEM) grid used to estimate the elevation of each raster cell must be downloaded manually from PRISM. Visit http://www.prism.oregonstate.edu/normals/, download the elevation grid in .BIL format at 4km and 800m spatial resolution.

# Transfer manually downloaded DEM files to working directory using interactive dialogue (if not previously transfered).
# This section will prompt the user to select the .ZIP file containing the 4km resolution DEM, then the .ZIP file containing the 800m resulution DEM.
if(!file.exists(paste0(dnld.dir, "Elevation/PRISM 4km DEM"))){
    print("Choose the 4km resolution DEM from the download location")
    dir.create(paste0(dnld.dir, "Elevation/PRISM 4km DEM"), recursive=T)
    unzip(file.choose(), exdir=paste0(dnld.dir, "Elevation/PRISM 4km DEM/"))
}
if(!file.exists(paste(dnld.dir, "Elevation/PRISM 800m DEM", sep="/"))){
    print("Choose the 800m resolution DEM from the download location")
    dir.create(paste0(dnld.dir, "Elevation/PRISM 800m DEM"), recursive=T)
    unzip(file.choose(), exdir=paste0(dnld.dir, "Elevation/PRISM 800m DEM/"))
}

## 4.2. Read in all rasters, crop to focal region, and calculate precip across months ====
PrecipSp.files <- as.list(gsub(".zip", ".bil", PrecipSp.files))
PrecipSp.rasters <- stack(lapply(PrecipSp.files, function(x){crop(raster(x), focal.region)}))
PrecipSp.mean <- overlay(PrecipSp.rasters, fun=mean)

if(!exists("elev800m.raster")){
    elev800m.raster <- crop(raster(paste0(dnld.dir, 
                                          "Elevation/PRISM 800m DEM/PRISM_us_dem_800m_bil.bil")), 
                            focal.region) 
}

## 4.3. Extract spatially paired values from precip and elevation raster ====
PrecipSp.df <- data.frame(PrecipSp.mean=PrecipSp.mean@data@values, elev=elev800m.raster@data@values)

## 4.4. Run regression analysis to test for change in precip over elevation ====
PrecipSp.mod <- lm(PrecipSp.mean~elev, data=PrecipSp.df)
summary(PrecipSp.mod)

## 4.5. Plot ====
gg.PrecipSp <- ggplot(PrecipSp.df, aes(elev, PrecipSp.mean))+
    geom_point(size=8, shape=16, color="grey50")+
    scale_x_continuous(name="Elevation (m)", limits=c(1800, 4250), breaks=c(2000, 3000, 4000), expand=c(0, 0))+
    scale_y_continuous(name="Total precip. (mm)", limits=c(20, 80), breaks=c(20, 40, 60, 80), 
                       expand=c(0, 0))+
    geom_smooth(method="lm", size=3, se=F, color="dodgerblue3")+
    theme_classic()
gg.PrecipSp

#########################################################
## 5. Precipitation over time ----
#########################################################
## Datasets: NOAA NCDC GHCN v2 Monthly Summaries
##           USDA NRCS WCC SNOTEL
## Variable: Total monthly precipitation (not available in GHCNv3 as of Feb 2015)
## URL: http://www.ncdc.noaa.gov/ghcnm/v2.php
##      http://www.wcc.nrcs.usda.gov/snow/

## 5.1. Download, import, clean data ====
# Create directory and download documentation
if(!file.exists(paste0(dnld.dir, "Precipitation/NOAA GHCN precip"))){
    dir.create(paste0(dnld.dir, "Precipitation/NOAA GHCN precip"), recursive=T)
    download.file("http://www1.ncdc.noaa.gov/pub/data/ghcn/v2/v2.prcp.readme", 
                  paste0(dnld.dir, "Precipitation/NOAA GHCN precip/v2.prcp.readme.txt"), quiet=T)
    download.file("http://www1.ncdc.noaa.gov/pub/data/cdo/documentation/GHCNDMS_documentation.pdf", 
                  paste0(dnld.dir, "Precipitation/NOAA GHCN precip/GHCNDMS_documentation.pdf"), quiet=T) 
}
if(!file.exists(paste0(dnld.dir, "Precipitation/USDA NRCS WCC SNOTEL"))){
    dir.create(paste0(dnld.dir, "Precipitation/USDA NRCS WCC SNOTEL"), recursive=T)
    download.file("http://www.wcc.nrcs.usda.gov/snotel/SNOTEL_brochure.pdf", 
                  paste0(dnld.dir, "Precipitation/USDA NRCS WCC SNOTEL/SNOTEL_brochure.pdf"), quiet=T)
}

# Download, import, and clean monthly data + metadata -- GHCN
if(!file.exists(paste0(dnld.dir, "Precipitation/NOAA GHCN precip/v2.prcp_adj"))){
    download.file("ftp://anonymous@ftp.ncdc.noaa.gov/pub/data/ghcn/v2/v2.prcp_adj.Z", 
                  paste0(dnld.dir, "Precipitation/NOAA GHCN precip/v2.prcp_adj.Z"), 
                  mode="w")
    download.file("ftp://ftp.ncdc.noaa.gov/pub/data/ghcn/v2/v2.prcp.inv", 
                  paste0(dnld.dir, "Precipitation/NOAA GHCN precip/v2.prcp.inv"), 
                  mode="w")
    setwd(paste0(dnld.dir, "Precipitation/NOAA GHCN precip/"))
    system2("uncompress", args="v2.prcp_adj.Z")
}

ghcn.Precip.df <- read.fwf(paste0(dnld.dir, "Precipitation/NOAA GHCN precip/v2.prcp_adj"), 
                           widths=c(11, 1, 4, rep(5, 12)), 
                           col.names=c("ID", "DUPLI", "YEAR", 
                                       paste0(rep("VALUE", times=12), 1:12)), 
                           header=F, stringsAsFactors=F)
if(!exists("ghcn.Precip.meta")){
    ghcn.Precip.meta <- read.fwf(file=paste0(dnld.dir, "Precipitation/NOAA GHCN precip/v2.prcp.inv"), 
                                 widths=c(11, 21, 11, 7, 8, 5), 
                                 col.names=c("ID", "NAME", "COUNTRY", "LATITUDE", "LONGITUDE", "STNELEV"), 
                                 header=F, stringsAsFactors=F, strip.white=T, comment.char="~")
}
ghcn.Precip.df[ghcn.Precip.df=="-9999"] <- NA # set missing data values

# Download, import, and clean monthly data + metadata -- SNOTEL
if(!file.exists(paste0(dnld.dir,"Precipitation/USDA NRCS WCC SNOTEL/SNOTEL_precip_", 
                       focal.years[1], "-", focal.years[2], ".csv"))){
    download.file(paste0("http://www.wcc.nrcs.usda.gov/reportGenerator/view_csv/customMultipleStationReport,metric/monthly/state=%22CO%22%20AND%20network=%22SNTLT%22,%22SNTL%22%7Cname/",
                         focal.years[1], "-01-01,", focal.years[2], 
                         "-12-31/stationId,elevation,latitude,longitude,PREC::value:monthly%20SUM"),
                  destfile=paste0(dnld.dir, "Precipitation/USDA NRCS WCC SNOTEL/SNOTEL_precip_", 
                                  focal.years[1], "-", focal.years[2], ".csv"), cacheOK=F)
}

snotel.Precip.df <- read.csv(paste0(dnld.dir, "Precipitation/USDA NRCS WCC SNOTEL/SNOTEL_precip_", 
                                    focal.years[1], "-", focal.years[2], ".csv"), skip=5, header=T, 
                             stringsAsFactors=F)
snotel.Precip.df <- with(snotel.Precip.df, data.frame(ID=Station.Id, 
                                                      YEAR=as.numeric(substring(Date, 5)), 
                                                      MONTH=sapply(substring(Date, 1, 3), 
                                                                   function(x)match(tolower(x), 
                                                                                    tolower(month.abb))), 
                                                      LATITUDE=Latitude, 
                                                      LONGITUDE=Longitude, 
                                                      ELEVATION=round(Elevation..ft./3.28084), 
                                                      MEAN.VALUE=Sum.Monthly.Precipitation.Accumulation..mm.))

## 5.2. Subset to stations within focal time period and focal region ====
# To focal time period
ghcn.Precip.df.crop <- subset(ghcn.Precip.df, subset=YEAR%in%focal.years[1]:focal.years[2], 
                              select=c("ID", "YEAR", paste0(rep("VALUE", times=length(focal.months)), 
                                                            focal.months)))
snotel.Precip.df.crop <- subset(snotel.Precip.df, subset=snotel.Precip.df$MONTH%in%focal.months)

# To focal region (exclude stations with coordinates outside of focal region)
ghcn.Precip.spdf <- with(ghcn.Precip.meta,
                         SpatialPointsDataFrame(coords=data.frame(LONGITUDE, LATITUDE), 
                                                data=ghcn.Precip.meta[, c(1:3, 5)]))
ghcn.Precip.df.focal <- ghcn.Precip.df.crop[ghcn.Precip.df.crop$ID%in%
                                                ghcn.Precip.spdf[!is.na(ghcn.Precip.spdf%over%
                                                                            as(focal.region*region.mult, 
                                                                               "SpatialPolygons")), ]@data[, 1], ]
ghcn.Precip.spdf.focal <- ghcn.Precip.spdf[ghcn.Precip.spdf$ID%in%
                                               ghcn.Precip.spdf[!is.na(ghcn.Precip.spdf%over%
                                                                           as(focal.region*region.mult, 
                                                                              "SpatialPolygons")), ]@data[, 1], ]
ghcn.Precip.df.focal$ID <- factor(ghcn.Precip.df.focal$ID)
# report included station metadata
ghcn.Precip.meta[ghcn.Precip.meta$ID%in%ghcn.Precip.df.focal$ID, 1:5]

snotel.Precip.spdf <- with(snotel.Precip.df.crop[!duplicated(snotel.Precip.df.crop$ID), ], 
                           SpatialPointsDataFrame(coords=data.frame(LONGITUDE, LATITUDE), 
                                                  data=snotel.Precip.df.crop[!duplicated(snotel.Precip.df.crop$ID), ]))
snotel.Precip.df.focal <- snotel.Precip.df.crop[snotel.Precip.df.crop$ID%in%
                                                    snotel.Precip.spdf[!is.na(snotel.Precip.spdf%over%
                                                                                  as(focal.region*region.mult, 
                                                                                     "SpatialPolygons")), ]@data$ID, ]
snotel.Precip.df.focal$ID <- factor(snotel.Precip.df.focal$ID)

## 5.3. Convert data to whole units (default is mm * 10) ====
cols <- paste0(rep("VALUE", times=length(focal.months)), focal.months)
ghcn.Precip.df.focal[cols] <- ghcn.Precip.df.focal[cols]/10
snotel.Precip.df.focal$MEAN.VALUE <- snotel.Precip.df.focal$MEAN.VALUE/10

## 5.4. Calculate mean of climate variable across focal months ====
ghcn.Precip.df.focal$MEAN.VALUE <- rowMeans(sapply(as.list(paste0("VALUE", 
                                                                  focal.months[1]:focal.months[2])), 
                                                   get, pos=ghcn.Precip.df.focal), 
                                            na.rm=F)
# Aggregate into season mean
snotel.Precip.df.focal <- aggregate(MEAN.VALUE~ID+YEAR+ELEVATION+LATITUDE+LONGITUDE, snotel.Precip.df.focal, 
                                    mean, na.rm=F, na.action=NULL)

## 5.5. Run regression analysis to test for change in precip over time ====
# COMBINED GHCN+SNOTEL DATA
PrecipTime.df <- merge(ghcn.Precip.df.focal, snotel.Precip.df.focal, all=T)
PrecipTime.df$NETWORK <- as.factor(c(rep("GHCN", times=nrow(ghcn.Precip.df.focal)), 
                                     rep("SNOTEL", times=nrow(snotel.Precip.df.focal))))
PrecipTime.lm <- lm(MEAN.VALUE~YEAR+ID, data=PrecipTime.df)
summary(PrecipTime.lm)
Anova(PrecipTime.lm, type="III")
partial.R2(update(PrecipTime.lm, ~. -YEAR), PrecipTime.lm)
PrecipTime.resids <- with(effect("YEAR", PrecipTime.lm, partial.residuals=T), 
                          data.frame(YEAR=data$YEAR, RESID=residuals, row.names=NULL))
PrecipTime.resids <- PrecipTime.resids %>%
    group_by(YEAR) %>%
    summarize(RESID=mean(RESID))

## 5.6. Plot ====
gg.PrecipTime <- ggplot(PrecipTime.resids, aes(x=YEAR, y=RESID))+
    geom_point(size=8, shape=16, color="grey50")+
    geom_smooth(method="lm", size=3, se=F, color="dodgerblue3")+
    scale_x_continuous(name="Year", limits=focal.years, breaks=c(1980, 1990, 2000, 2010), expand=c(0, 0))+
    scale_y_continuous(name="Total precipitation (mm)", limits=c(20, 80), breaks=c(20, 40, 60, 80),
                       expand=c(0, 0))+
    theme_classic()
gg.PrecipTime

#########################################################
## 6. Snowmelt over space ----
#########################################################
## Dataset: Personal HOBO loggers + SNOTEL
## Variable: DOY snowmelt (measured as first day where snow depth is 0 cm)
## URL: http://www.wcc.nrcs.usda.gov/snow/

## 6.1. Import data ====
# Looks for file 'SnowMeltSpace.csv' from Dryad repository in 'dnld.dir' specified in section 1.2
SnwSp.data <- read.csv(paste0(dnld.dir, "SnowMeltSpace.csv"), header=T)

## 6.2. Run regression analysis to test for change in snow melt over elevation ====
SnwSp.mod <- lm(melt.doy~elevation, data=SnwSp.data)
summary(SnwSp.mod)

## 6.3. Plot ====
gg.SnwSp <- ggplot(SnwSp.data, aes(x=elevation, y=melt.doy))+
    geom_point(size=8, shape=16, color="grey50")+
    geom_smooth(data=SnwSp.data, aes(x=elevation, y=melt.doy), method="lm", size=3, 
                se=F, color="dodgerblue3")+
    scale_x_continuous(name="Elevation (m)", limits=c(1800, 4250), breaks=c(2000, 3000, 4000),
                       expand=c(0, 0))+
    scale_y_continuous(name="Snowmelt DOY", limits=c(90, 180), breaks=seq(90, 180, 30),
                       expand=c(0, 0))+
    theme_classic()
gg.SnwSp

#########################################################
## 7. Snowmelt over time ----
#########################################################
## Dataset: billy barr Snowmelt Dataset
## Variable: DOY snowmelt
## URL: http://www.gothicwx.com/

## 7.1. Import data ====
# Looks for file 'SnowMeltTime.csv' from Dryad repository in 'dnld.dir' specified in section 1.2
SnwTime.data <- read.csv(paste0(dnld.dir, "SnowMeltTime.csv"), header=T)

## 7.2. Run regression analysis to test for change in snow melt over time ====
SnwTime.mod <- lm(gothic.melt.doy~year, data=SnwTime.data)
summary(SnwTime.mod)
summary(effect("year", SnwTime.mod))
# average change over time
max(summary(effect("year", SnwTime.mod, xlevels=37))$effect)-
    min(summary(effect("year", SnwTime.mod, xlevels=37))$effect)

## 7.3. Plot ====
gg.SnwTime <- ggplot(SnwTime.data, aes(x=year, y=gothic.melt.doy))+
    geom_point(size=8, shape=16, color="grey50")+
    geom_smooth(data=SnwTime.data, aes(x=year, y=gothic.melt.doy), method="lm", size=3, 
                se=F, linetype="dashed", color="dodgerblue3")+
    scale_x_continuous(name="Year", limits=focal.years, breaks=c(1980, 1990, 2000, 2010),
                       expand=c(0, 0))+
    scale_y_continuous(name="Snowmelt DOY", limits=c(90, 180), breaks=seq(90, 180, 30),
                       expand=c(0, 0))+
    theme_classic()
gg.SnwTime

#########################################################
## 8. Soil moisture over space ----
#########################################################
## Dataset: USDA NRCS SNOTEL + RMBLNet
## Variable: Mean monthly soil moisture
## URL: http://www.wcc.nrcs.usda.gov/snow/
##      http://www.wrcc.dri.edu/rmbl/

## 8.1. Download, import, clean data ====
# Specify stations
# 1°x1° stations include RMBLNet + SNOTEL 369, 1141, 680, 380, 737
# 2°x2° stations include 1°x1° stations + SNOTEL 762, 1059, 701, 682, 378, 939, 531, 938, 485, 415, 802, 1014

snotel.stations <- c(369, 1141, 680, 380, 737, 762, 1059, 701, 682, 378, 939, 531, 938, 485, 415,
                     802, 1014)
snotel.elevs <- data.frame(Station.Id=snotel.stations, Elevation=c(3231L, 3243L, 2926L, 3097L, 3261L,
                                                                   3487L, 3054L, 3280L, 3036L, 2865L,
                                                                   3158L, 3474L, 3399L, 3474L, 3216L,
                                                                   2865L, 2725L))
SoilSp.url <- paste0("http://www.wcc.nrcs.usda.gov/reportGenerator/view_csv/customMultipleStationReport,metric/monthly/", paste0(snotel.stations, collapse=":CO:SNTL%7C"), ":CO:SNTL%7Cid=%22%22%7Cname/2009-06-01,2015-09-30/stationId,name,SMS:-2:value:monthly%20MEAN,SMS:-4:value:monthly%20MEAN,SMS:-8:value:monthly%20MEAN,SMS:-20:value:monthly%20MEAN,SMS:-40:value:monthly%20MEAN")

# Create directory and download documentation
if(!file.exists(paste0(dnld.dir, "Soil Moisture/USDA NRCS WCC SNOTEL"))){
    dir.create(paste0(dnld.dir, "Soil Moisture/USDA NRCS WCC SNOTEL"), recursive=T)
    download.file("http://www.wcc.nrcs.usda.gov/snotel/SNOTEL_brochure.pdf", 
                  paste0(dnld.dir, "Soil Moisture/USDA NRCS WCC SNOTEL/SNOTEL_brochure.pdf"), quiet=T)
}

# Download and import monthly SNOTEL soil moisture data
if(!file.exists(paste0(dnld.dir, "Soil Moisture/USDA NRCS WCC SNOTEL/soilmoisture.csv"))){
    download.file(SoilSp.url, paste0(dnld.dir, "Soil Moisture/USDA NRCS WCC SNOTEL/soilmoisture.csv"))
}

# Import SNOTEL DATA + clean
SnwSp.snotel <- tbl_df(read.csv(paste0(dnld.dir, 
                                       "Soil Moisture/USDA NRCS WCC SNOTEL/soilmoisture.csv"), skip=5))
SnwSp.snotel <- SnwSp.snotel %>%
    rename(soilm5cm=Mean.Monthly.Soil.Moisture.Percent..2in..pct., 
           soilm10cm=Mean.Monthly.Soil.Moisture.Percent..4in..pct., 
           soilm20cm=Mean.Monthly.Soil.Moisture.Percent..8in..pct., 
           soilm50cm=Mean.Monthly.Soil.Moisture.Percent..20in..pct., 
           soilm100cm=Mean.Monthly.Soil.Moisture.Percent..40in..pct.) %>%
    mutate(soilm5cm=ifelse(soilm5cm<100, soilm5cm/100, NA), 
           soilm10cm=ifelse(soilm10cm<100, soilm10cm/100, NA), 
           soilm20cm=ifelse(soilm20cm<100, soilm20cm/100, NA), 
           soilm50cm=ifelse(soilm50cm<100, soilm50cm/100, NA), 
           soilm100cm=ifelse(soilm100cm<100, soilm100cm/100, NA)) %>%
    separate("Date", c("Month", "Year"), sep=" ", convert=T) %>%
    mutate(Month=match(Month, month.abb), dsource="SNOTEL")
SnwSp.snotel <- full_join(SnwSp.snotel, snotel.elevs, by="Station.Id") %>%
    dplyr::select(-Station.Id) # specify package to avoid conflict with raster::select

# Read RMBLNet Data + clean
SnwSp.rmbl <- tbl_df(read.csv(paste0(dnld.dir, "RMBLNetSoilMoisture.csv")))
SnwSp.rmbl <- SnwSp.rmbl %>%
    dplyr::select(-contains("min")) %>%
    dplyr::select(-contains("max")) %>%
    rename(soilm5cm=SoilM5cm.avg, soilm25cm=SoilM25cm.avg, soilm50cm=SoilM50cm.avg) %>%
    mutate(soilm5cm=ifelse(soilm5cm<1, soilm5cm, NA), 
           soilm25cm=ifelse(soilm25cm<1, soilm25cm, NA), 
           soilm50cm=ifelse(soilm50cm<1, soilm50cm, NA)) %>%
    separate("Date", c("Month", "Day", "drop"), sep="/", convert=T) %>%
    dplyr::select(-drop, -DOY, -DOR) %>%
    group_by(Station.Name, Elevation, Month, Year) %>%
    summarize(soilm5cm=mean(soilm5cm, na.rm=T), soilm25cm=mean(soilm25cm, na.rm=T), 
              soilm50cm=mean(soilm50cm, na.rm=T)) %>%
    mutate(dsource="RMBL")

# Merge data sources (warns about joining factors with different levels)
SnwSp.full <- full_join(SnwSp.snotel, SnwSp.rmbl) %>%
    mutate(Station.Name=as.factor(Station.Name), dsource=as.factor(dsource)) %>%
    mutate(soilm2025cm=pmax(soilm20cm, soilm25cm, na.rm=T))

## 8.2. Subset to stations within focal time period & average across growing season ====
SnwSp.full %<>%
    filter(Month %in% 6:9) %>%
    group_by(Station.Name, Elevation, dsource) %>%
    summarize(soilm2025cm=mean(soilm2025cm, na.rm=T)) %>%
    filter(Station.Name!="Middle Fork Camp") # remove replicate where deeper soil probes aren't working

# 8.3. Run regression analysis to test for change in soil moisture over elevation ====
SnwSp.mod <- lm(soilm2025cm~Elevation, data=SnwSp.full)
summary(SnwSp.mod)
plot(effect("Elevation", SnwSp.mod, partial.residuals=T), partial.residuals="raw", smooth.residuals=F)

# 8.4. Plot ====
gg.SoilSp <- ggplot(SnwSp.full, aes(x=Elevation, y=soilm2025cm))+
    geom_point(size=8, shape=16, color="grey50")+
    geom_smooth(size=3, method="lm", se=F, linetype="dashed", color="dodgerblue3")+
    scale_x_continuous(name="Elevation (m)", limits=c(1800, 4250), breaks=c(2000, 3000, 4000), expand=c(0, 0)) +
    scale_y_continuous(name="Mean soil moist. (vol%)", limits=c(0, 0.5), breaks=c(0, 0.25, 0.5), labels=c(0, 25, 50), expand=c(0, 0))+
    theme_classic()
gg.SoilSp

#########################################################
## 9. Soil moisture over time ----
#########################################################
## Datasets: NOAA Climate Prediction Center soil moisture
## Variable: Mean monthly soil moisture
## URL: ftp://ftp.cdc.noaa.gov/Datasets/cpcsoil/

## 9.1 Download, import, clean data ====
# Create directory and download documentation
if(!file.exists(paste0(dnld.dir, "Soil Moisture/NOAA CPC soilm"))){
    dir.create(paste0(dnld.dir, "Soil Moisture/NOAA CPC soilm"), recursive=T)
    download.file("ftp://ftp.cdc.noaa.gov/Datasets/cpcsoil/README", 
                  paste0(dnld.dir, "Soil Moisture/NOAA CPC soilm/README.txt"), quiet=T)
    download.file("http://www.cpc.noaa.gov/products/Soilmst_Monitoring/Papers/2003JD004345.pdf", 
                  paste0(dnld.dir, "Soil Moisture/NOAA CPC soilm/Fan&2004_JGR.pdf"), quiet=T)
    download.file("http://www.cpc.ncep.noaa.gov/products/Soilmst_Monitoring/Papers/jhuang.pdf", 
                  paste0(dnld.dir, "Soil Moisture/NOAA CPC soilm/Huang&1996_JC.pdf"), quiet=T)
}

# Download and import monthly data
if(!file.exists(paste0(dnld.dir, "Soil Moisture/NOAA CPC soilm/soilw.mon.mean.v2.nc"))){
    download.file("ftp://ftp.cdc.noaa.gov/Datasets/cpcsoil/soilw.mon.mean.v2.nc", 
                  paste0(dnld.dir, "Soil Moisture/NOAA CPC soilm/soilw.mon.mean.v2.nc"), quiet=T)
}

SoilTime.stack <- stack(paste0(dnld.dir, "Soil Moisture/NOAA CPC soilm/soilw.mon.mean.v2.nc"))

## 9.2. Subset to stations within focal time period and focal region ====
# To focal time period
SoilTime.sub <- which(as.numeric(substr(names(SoilTime.stack), 2, 5))%in%(focal.years[1]:focal.years[2])&
                          as.numeric(substr(names(SoilTime.stack), 7, 8))%in%focal.months)
SoilTime.stack <- subset(SoilTime.stack, SoilTime.sub)

# To focal region (first convert focal region to 360 coordinate reference system)
focal.region360 <- focal.region
focal.region360@xmin <- 360+focal.region360@xmin
focal.region360@xmax <- 360+focal.region360@xmax
# Exclude stations with coordinates outside of focal region
SoilTime.sub.stack <- crop(SoilTime.stack, focal.region360)

# 9.3. Convert data to volumetric percentage ====
# following methods of http://www.fao.org/docrep/r4082e/r4082e03.htm using 760mm soil column per NOAA CPC documentation
SoilTime.sub.stack <- stack((SoilTime.sub.stack/1000))

# 9.4. Calculate mean soil moisture across growing season ====
SoilTime.means <- colSums(raster::extract(SoilTime.sub.stack, 1:4))/4 # preface with package name to avoid conflict with tidyr::extract
SoilTime.data <- data.frame(year=as.numeric(unique(substr(names(SoilTime.means), 2, 5))), 
                            year.means=sapply(split(SoilTime.means, 
                                                    (seq_along(SoilTime.means)-1)%/%length(focal.months)),
                                              mean))

# 9.5. Run regression analysis to test for change in soil moisture over time ====
SoilTime.mod <- lm(year.means~year, data=SoilTime.data)
summary(SoilTime.mod)

# 9.6. Plot ====
gg.SoilTime <- ggplot(SoilTime.data, aes(x=year, y=year.means))+
    geom_point(size=8, shape=16, color="grey50")+
    geom_smooth(size=3, method="lm", se=F, color="dodgerblue3")+
    scale_x_continuous(name="Year", limits=focal.years, breaks=c(1980, 1990, 2000, 2010),
                       expand=c(0, 0))+
    scale_y_continuous(name="Mean soil moist. (vol%)", limits=c(0, 0.5), breaks=c(0, 0.25, 0.5), 
                       labels=c(0, 25, 50), expand=c(0, 0))+
    theme_classic()
gg.SoilTime

#########################################################
## 10. Calculate paces of climate change ----
#########################################################
## This section estimates the pace of climate change for temperature, precipitation, snow melt, and soil moisture. The endpoints of each variable's regression over time are compared to the regression over space to estimate the elevational shift of that climate indicator. Dividing by the intervening number of decades produces a pace of climate change (m/decade). The standard errors of each regression are propogated to robustly estimate the uncertainty around the pace of change in each variable.

## 10.1. Temperature pace ====
TmeanSp.summary <- summary(TmeanSp.mod) # examine regression of temperature over space
TmeanTime.summary <- effect("YEAR", ghcn.summary.lm, xlevels=list(YEAR=c(1978, 2014)))
Tmean.params <- data.frame(sp_int=TmeanSp.summary$coefficients["(Intercept)", 1:2], 
                           sp_b=TmeanSp.summary$coefficients["elev", 1:2], 
                           contemp=c(TmeanTime.summary$fit[2], TmeanTime.summary$se[2]), 
                           histor=c(TmeanTime.summary$fit[1], TmeanTime.summary$se[1]))
Tmean.expr <- expression(((histor-sp_int)/sp_b)-((contemp-sp_int)/sp_b))
Tmean.res <- propagate(Tmean.expr, Tmean.params, nsim=1e6)
mean(Tmean.res$resSIM/(2014-1978))*10 # mean
(Tmean.res$sim[c("2.5%", "97.5%")]/(2014-1978))*10 # lower/upper 95% CIs
(diff(quantile(Tmean.res$resSIM, c(0.025, 0.975)))/(201.4-197.8))/(qnorm(0.975)*2) # standard error in units of m/decade

## 10.2. Precipitation pace ====
PrecipSp.summary <- summary(PrecipSp.mod)
PrecipTime.summary <- effect("YEAR", PrecipTime.lm, xlevels=list(YEAR=c(1978, 2014)))

Precip.params <- data.frame(sp_int=PrecipSp.summary$coefficients["(Intercept)", 1:2], 
                            sp_b=PrecipSp.summary$coefficients["elev", 1:2], 
                            contemp=c(PrecipTime.summary$fit[2], PrecipTime.summary$se[2]), 
                            histor=c(PrecipTime.summary$fit[1], PrecipTime.summary$se[1]))
Precip.expr <- expression(((histor-sp_int)/sp_b)-((contemp-sp_int)/sp_b))
Precip.res <- propagate(Precip.expr, Precip.params, nsim=1e6)
mean(Precip.res$resSIM/(2014-1978))*10 # mean
(Precip.res$sim[c("2.5%", "97.5%")]/(2014-1978))*10 # lower/upper 95% CIs
(diff(quantile(Precip.res$resSIM, c(0.025, 0.975)))/(201.4-197.8))/(qnorm(0.975)*2) # standard error in units of m/decade

## 10.3. Snowmelt pace ====
SnwSp.summary <- summary(SnwSp.mod)
SnwTime.summary <- effect("year", SnwTime.mod, xlevels=list(year=c(1978, 2014)))

Snw.params <- data.frame(sp_int=SnwSp.summary$coefficients["(Intercept)", 1:2], 
                         sp_b=SnwSp.summary$coefficients["elevation", 1:2], 
                         contemp=c(SnwTime.summary$fit[2], SnwTime.summary$se[2]), 
                         histor=c(SnwTime.summary$fit[1], SnwTime.summary$se[1]))
Snw.expr <- expression(((histor-sp_int)/sp_b)-((contemp-sp_int)/sp_b))
Snw.res <- propagate(Snw.expr, Snw.params, nsim=1e6)
mean(Snw.res$resSIM/(2014-1978))*10 # mean
(Snw.res$sim[c("2.5%", "97.5%")]/(2014-1978))*10 # lower/upper 95% CIs
(diff(quantile(Snw.res$resSIM, c(0.025, 0.975)))/(201.4-197.8))/(qnorm(0.975)*2) # standard error in units of m/decade

## 10.4. Soil moisture pace ====
SoilSp.summary <- summary(SoilSp.mod)
SoilTime.summary <- effect("year", SoilTime.mod, xlevels=list(year=c(1978, 2014)))

Soil.params <- data.frame(sp_int=Soil.space.summary$coefficients["(Intercept)", 1:2], 
                          sp_b=Soil.space.summary$coefficients["Elevation", 1:2], 
                          contemp=c(SoilTime.summary$fit[2], SoilTime.summary$se[2]), 
                          histor=c(SoilTime.summary$fit[1], SoilTime.summary$se[1]))
Soil.expr <- expression(((histor-sp_int)/sp_b)-((contemp-sp_int)/sp_b))
Soil.res <- propagate(Soil.expr, Soil.params, nsim=1e6)
mean(Soil.res$resSIM/(2014-1978))*10 # mean
(diff(quantile(Soil.res$resSIM, c(0.025, 0.975)))/(201.4-197.8))/(qnorm(0.975)*2) # standard error in units of m/decade