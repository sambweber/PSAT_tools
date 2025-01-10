#####################################################################################################
#CONTAINS A SET OF CUSTOM FUNCTIONS FOR WORKING WITH OUTPUTS FILES FROM WILDLIFE COMPUTERS MINIPAT TAGS
#####################################################################################################

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# code by Sam Weber, Ascension Island Government Conservation & Fisheries Department, June 2016
# sam.weber@exeter.ac.uk
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' Contains functions for processing NetCDF files output from the Wildlife Computers GPE3 model
#' and extracting volume contours describing the residency areas of individuals and populations

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# volume_contour
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' NOTE: Elements of this function were adapted from code contibuted by Josh Stewart from Scripps Institue of Oceanography
#' 
#' Extracts precentage volume contours from GPE3 likelihood rasters at the specified \code{levels}. 
#' These can be interpreted as the smallest area in which the animal had a given probability of being located during the time
#' period represented in the \code{input} layer. The function is primarily intended to be used within the \code{gpe_residency}
#' script below.

#' @param input a \code{Raster} object containing gridded likelihoods output by GPE3 (either a single
#     twelve hourly surface or a composited layer for an individual animal or population)
#' @param res.out a scaling factor used to resample the original likelihood raster in order to generate smooth contours. 
#' Before contouring, the \code{input} Raster is bilinearly resampled onto a grid at \code{res.out} times the original resolution.
#' @param levels a numeric vector of percentages corresponding to the volume contours required.
#' @return Returns a three element list containing the resampled \code{Raster} and a \code{SpatialLinesDataFrame} and \code{SpatialPolygonsDataFrame}
#' for the volume contours.

volume_contour <- function(input, res.out = 10, levels = c(95,75,25,50),reclassify = T){
  
  require(raster); require(rgdal); require(sf); require(rgeos)
  
  resamp = raster::disaggregate(x=input,fact = res.out,method="bilinear")
  Z = raster::values(resamp)
  resamp[] = Z/sum(Z,na.rm=T) # normalise 
  
  vclevels = sort(levels,decreasing = T) #sort the volume contour levels
  rastvals = sort(raster::values(resamp)) #rank the raster probability values
  breaks   = sapply(1-vclevels/100,FUN=function(x){rastvals[max(which(cumsum(rastvals)<=x))]})
  breaks   = c(breaks,1)
  
  cols = rev(heat.colors(length(levels)))
  plot(resamp,col=cols,breaks=breaks,legend = F,xlab="Longitude",ylab="Latitude",bg=400)
  legend("topright",legend = paste(vclevels,"%",sep=" "),fill = cols)
  
  #add contour lines
  vcs = rasterToContour(resamp,levels = breaks[1:length(levels)],maxpixels = ncell(resamp))
  vcs$level = levels
  plot(vcs,add=T)
  
  #convert to polygons - at the moment conversion to 'Spatial' class fails when there are multiple polygons in a level using st_polygonize 
  #or when there is only one polygon in a level using st_cast so adopting a 'try' approach instead until this issue is resolved
  #vcpolys = try(st_as_sf(vcs) %>% st_cast('POLYGON') %>% as('Spatial'),silent=TRUE) 
  #if(class(vcpolys) == "try-error") vcpolys = st_as_sf(vcs) %>% st_polygonize() %>% as('Spatial')
  #vcpolys = rgeos::createSPComment(vcpolys)
  vcpolys = st_as_sf(vcs) %>% st_cast('POLYGON')
  
  #calculate area using Lambert Azimuthal Equal Area projection
  centroid = st_coordinates(st_centroid(vcpolys))
  lamberts = paste0("+proj=laea +lon_0=",centroid[1,1]," +lat_0=",centroid[1,2])
  vcpolys$area_km2 = as.numeric(round(st_area(st_transform(vcpolys,lamberts))/1000000,1))
  
  #reclassify raster by levels
  if(reclassify){
  resamp = reclassify(resamp,matrix(c(min(rastvals),breaks[-length(breaks)],breaks,NA,vclevels),ncol=3),include.lowest=T)
  }
  
  #return outputs as a list
  return(list(resamp,vcs,vcpolys))
  
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# gpe_residency
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' A function for calculating the mean residency likelihoods and associated volume contours
#' for an individual or population across a deployment period. This function is a wrapper
#' for the volume_contour script above. It first calculates the mean 12 hourly residency 
#' likelihoods for the individual or population on a standard grid (using an optional weighting
#' for track length) and then extracts percentage volume contours at the levels specified 

#' @param gpe3_files a character vector containing file paths to the GPE3.nc files for a collection of animals.
#'  If \code{gpe3_files} contains a single element, the residency surface for a single animal will be returned.
#'  The \code{weight} argument is meaningless in this case. 
#'  
#' @param weight a logical argument indicating whether or not to use a weighted average for calculating
#'  mean residency likelihoods. If \code{weight = TRUE} each 12 hourly-likelihood surface is multiplied 
#'  by the inverse of the number of individuals that had a measurement at this relative timepoint in their
#'  respective tracks. This effectively implements the method of Block et al.2011 and ensures that population 
#'  residency is not excessively influenced by large numbers of short term deployments close to the release site.
#'  
#'  @param percentile a numeric between 0 and 1 specifying the percentile of track duration at which weighting should
#'  be stopped. All 12-hourly surfaces for timepoints greater than this will be weighted equally.
#'  
#'  @param track_end an optional vector of class \code{POSIXct} containing the end points of each track
#'  which will be used to truncate the 12 hourly likelihood series. This can be useful for removing
#'  data for dead and hence stationary animals gathered before tags release (e.g. floaters) that result in artefacts 
#'  in residency surfaces. Vector should be equal in length to \code{gpe3_files}. If nothing is supplied then
#'  all 12 hourly surfaces will be used (the default)

#' @param ... Other arguments passed to \code{volume_contour}
#' 
#' @examples:
#  ncfiles = list.files('mydirectory',pattern='\\.nc$',full.names=T,recursive=T)
#  result  = gpe_residency(ncfiles,weight=T,levels=c(95,75,50,25))

gpe_residency <- function(gpe3_files, track_end = NULL, weight=F,percentile=0.85, ...){
  
  require(ncdf4)
  
  Z = lapply(gpe3_files,function(y) ncvar_get(nc_open(y),'twelve_hour_timestamps'))
  
  if(!is.null(track_end)){ 
    if (length(track_end) < length(Z)) {stop ("number of track end dates don't match number of files") 
  } else {Z = lapply(seq(Z),function(i) {which(Z[[i]] <= as.numeric(track_end[i]))})}
  } else {Z = lapply(Z,seq)}
  
  
  if (weight == T){
    
    counts = sapply(Z,length)  
    q = ceiling(quantile(counts,percentile))
    nbyday = sapply(seq(max(counts)),FUN = function(x){length(counts[counts>=x])})
    weights = c(1/nbyday[1:q],rep(1,length(nbyday) - q))
    
  }
  
  output = list()
  
  for(i in seq(gpe3_files)){
    
    likelihoods = raster::stack(gpe3_files[i],varname = "twelve_hour_likelihoods",bands=Z[[i]])
    if (weight == T) likelihoods = likelihoods*(weights[1:dim(likelihoods)[3]])
    av_likelihoods = overlay(likelihoods,fun=mean)
    av_likelihoods[is.na(av_likelihoods)] = 0
    output[[i]] = av_likelihoods
    
  }
  
  boundaries <- extent(c(min(do.call("cbind",lapply(output,FUN = function(x){extent(x)@xmin}))),
                         max(do.call("cbind",lapply(output,FUN = function(x){extent(x)@xmax}))),
                         min(do.call("cbind",lapply(output,FUN = function(x){extent(x)@ymin}))),
                         max(do.call("cbind",lapply(output,FUN = function(x){extent(x)@ymax})))))
  
  output = lapply(output,FUN = function(x){extend(x,boundaries,value=0)})
  
  population_average = overlay(do.call("stack",output),fun=mean) 
  
  return(volume_contour(population_average, ...))
  
  
}


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# gpe3_track_sf
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# returns trackpoints from GPE3 output as an sf object

# input is a directory containing GPE3 files. It assumes that deployment and release
# metadata have been entered on portal.

# If model_only = TRUE only the 12 hourly modelled locations from GPE3 are returned (not SST and light positions)
# If combine = TRUE, data from all individuals are combined into a single sf table

gpe3_track_sf = function(dirs,model_only=FALSE,combine=FALSE,include_endpoints = TRUE){
  
  result = lapply(dirs,function(dir){
  
  gpe3 = list.files(dir,full.names=TRUE,pattern = "\\GPE3.csv$") 
  
  if (length(gpe3) == 1L){
    
    gpe3 = read.csv(gpe3,skip=5,header=T,stringsAsFactors=FALSE) %>%
           dplyr::select(2:6) %>% 
           setNames(c("ptt","datetime","latitude","longitude","type")) %>%
           mutate(datetime = as.POSIXct(datetime,format = "%d-%b-%Y %H:%M:%S",tz="GMT")) %>%
           arrange(datetime) %>%
           mutate(type = gsub("Light -","",type)) %>%
           mutate(type = tolower(type)) %>%
           {if(model_only) dplyr::filter(.,type %in% c('none','user')) else .} %>%
           st_as_sf(coords = c('longitude','latitude'),crs=4326)
    
    usr = which(gpe3$type == "user")
    gpe3$type[usr[1L]] = "deployment"; gpe3$type[usr[length(usr)]] = "release"
    
    if(include_endpoints) return(gpe3) else return(subset(gpe3,!type %in% c('deployment','release')))
    
  } else {
    
    return(NULL)
    
  }
  
})
  
  setNames(result,basename(dirs)) %>% .[!vapply(., is.null, logical(1))] %>% {if(combine) do.call('rbind',.) else .}
  
}


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# gpe3_sim_track
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# This function simulates sets of N tracks by sampling from the 12-hourly 
# probability distributions output from GP3. Tracks can optionally be
# cut by a track_end (POSIXct) vector that is the same lenght as the GPE3
# files.

gpe3_sim_track = function(gpe3_files,track_end,N=100){
  
  source('SOFTWARE/R/custom_scripts/GridTools.R')
  
  if (length(track_end) != length(gpe3_files)) {stop ("number of track end dates don't match number of files")} 
  
  output = list()
  
  for(i in seq(gpe3_files)){
    
    Z =  ncvar_get(nc_open(gpe3_files[i]),'twelve_hour_timestamps') %>%
      as.POSIXct(origin='1970-01-01',tz='GMT')
    
    if(!is.null(track_end)){ 
      idx = which(Z <= track_end[i])
    } else {idx = lapply(Z,seq)}
    
    likelihoods = raster::stack(gpe3_files[i],varname = "twelve_hour_likelihoods",bands=idx)
    resamples   = sampleWeighted(likelihoods,N=N,sf=T) %>% 
      lapply(function(j) mutate(j,sim=1:N))
    output[[i]] = lapply(seq(resamples),function(j) mutate(resamples[[j]],index=j,datetime=Z[j])) %>%
      do.call(rbind,.)
    
  }
  
  return(output)
  
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# gpe3_meta
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Fetch metadata from a vector of GPE3 files
                         
gpe3_meta = function(files){

lapply(files,function(x){  
  
read.csv(x,nrows=1,skip=1,header=F,stringsAsFactors = FALSE)[1,] %>%  
  strsplit('@') %>% unlist() %>% 
  strsplit(',') %>% unlist() %>%
  .[-1]

})  
  
}


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# depth_profile
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Plots the depth series of a tag, optionally limiting to only the last 'n_tail' days to 
# help with diagnosing points of tag release, morality etc.

depth_profile = function(dir,n_tail){
  
  require(dplyr)
  
  series = list.files(dir,full.names=TRUE,pattern = "\\Series.csv$") %>% read.csv()
  summary = list.files(dir,full.names=TRUE,pattern = "\\Summary.csv$") %>% read.csv() %>% dplyr::filter(!is.na(Ptt))

  series = mutate(series,datetime = lubridate::dmy_hms(paste(Day,Time),tz='GMT'))
  if(!missing(n_tail)) series = subset(series,difftime(last(datetime),datetime,unit='days')<n_tail)

  plot(series$datetime,-(series$Depth),type='l',main=summary$Ptt)

}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# depth_profile
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# A function for plotting depth and temperature series data - the function allows user
# to interactively select a directory and then plots all series files within it in sequence

plot_series = function(){
dir = choose.dir()
f = list.files(dir,pattern = '*-Series.csv$',full.names = TRUE,recursive = TRUE) 

map(f, function(.x){
series = read.csv(.x) 
title = unique(series$Ptt) 

series = mutate(series, datetime = lubridate::dmy_hms(paste(Day,Time))) %>%
         dplyr::select(datetime,Temperature,Depth) %>%
         mutate(Depth = Depth * -1) %>%
         pivot_longer(cols = c(Temperature,Depth), names_to = 'variable')
 
ggplot(series,aes(x = datetime, y = value)) +
  geom_point(aes(colour = variable)) +
  facet_wrap(~variable,ncol=1,scales = 'free_y') + 
  labs(title = title)
})

}
                         

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# release
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Finds the location and time of the first argos transmission after release. Need to supply a directory 
# containing the files, or use the wcDownload to put these in a temporary folder.

# dir is the directory containing the WC portal files
# landmask is a spatialPolygon object e.g. obtained by gUnaryUnion(getMap("high")) from package rworldxtra

# If the drift_correct argument is specified then mean tag drift speed and bearing over 
# a period of 'drift_correct' hours after the first high quality location will be used to 
# back extrapolate the true location at the point the track ended. Note this is only 
# suitable when lags are fairly short (a few days maximum) and due to tags floating at
# the surface. For 'sitters' it may be more appropriate to use the lag from true release
# rather than the track end as animals sitting on the seabed are unlikely to drift at the
# same rate. For very long lags which may indicate the tag has been eaten(!) then the final
# geolocation before the track ended is probably the best bet.

release = function(dirs,landmask=NULL,interactive=F,truncate_hrs,n_tail=50,drift_correction){

require(sf); require(trip); require(lubridate)

result=list()  
  
for(d in seq(dirs)){  
    
summary = list.files(dirs[d],full.names=TRUE,pattern = "\\Summary.csv$") %>% read.csv() %>% dplyr::filter(!is.na(Ptt)) #somtimes summary table containsfirst row with no deployment data
argos = list.files(dirs[d],full.names=TRUE,pattern = "\\-Argos.csv$") %>% read.csv()

if(dim(summary)[1]>1) stop (paste('summary.csv for',dirs[d],'has more than 1 row!'))

if(is.na(summary$ReleaseDate)) {result[[d]] = data.frame("ptt" = summary$Ptt,
                                                         "datetime" = as.POSIXct(NA,tz='GMT'), 
                                                         "track_end" = as.POSIXct(NA,tz='GMT'),
                                                         "first_argos" = as.POSIXct(NA,tz='GMT'),
                                                         "longitude" = NA, 
                                                         "latitude" = NA,
                                                         "lag_release_hrs" = NA,
                                                         "lag_track_end_hrs" = NA,
                                                         "release_type" = 'non-transmitter')} else {
  
argos = argos[!is.na(argos$Latitude)&!is.na(argos$Longitude),]
argos$Date = parse_date_time(argos$Date,orders = c("HMS dbY","dmY HM"),tz='GMT')


# Filter out argos locations that are on land, low quality or have unrealistic drift speeds
#coordinates(argos) = ~Longitude + Latitude
#proj4string(argos) = CRS("+proj=longlat")
argos = st_as_sf(argos,coords = c('Longitude','Latitude'), crs=4326)

if(!is.null(landmask)){
  land = st_make_valid(land)
  argos <- argos[!st_intersects(argos,land,sparse=F)[,1],]
}
  
argos <- trip::forceCompliance(as(argos,'Spatial'),c("Date","Ptt"))
argos <- argos[speedfilter(argos,max.speed=3.5),]
argos <- as(argos,'SpatialPointsDataFrame') %>% st_as_sf()
argos <- argos[argos$LocationQuality %in% c("A","1","2","3"),]

# Find release time
release = parse_date_time(summary$ReleaseDate,orders = c("H:M:S d-b-Y","d/m/Y H:M"),tz='GMT')

# Argos locations after release time
argos <- argos[argos$Date > release,]

# Truncate the track by some automatic amount if the release was 
# premature

track_end = ifelse(

summary$ReleaseType %in% c('Premature','Floater','Sitter') & !missing(truncate_hrs),

release - 3600 * truncate_hrs, release

) %>% as.POSIXct(origin='1970-01-01',tz='GMT')
  

# If interactive
if(interactive){
  
  tryCatch(
    {depth_profile(dirs[d],n_tail = n_tail)
     abline(v=release,col='blue')
     abline(v=track_end,col='red')
     
     happy <- readline("Accept track end (Y or N)?") 
    
     while(happy == 'N'){
       depth_profile(dirs[d],n_tail = n_tail)
       abline(v=release,col='blue')
       print('Select the track end from the series')
       track_end = locator(1)$x
       abline(v=track_end,col='red')
       happy <- readline("Accept track end (Y or N)?") 
     }
    
     track_end = as.POSIXct(track_end,origin='1970-01-01',tz='GMT')  
    },error=function(e) print('No time series - using original track end')
  )
}  


t = difftime(argos$Date[1],release,units="hours")
t2 = difftime(argos$Date[1],track_end,units="hours")

result[[d]] = data.frame("ptt" = summary$Ptt,
                         "datetime" = release, 
                         "track_end" = track_end,
                         "first_argos" = argos$Date[1],
                         "longitude" = st_coordinates(argos)[1,1], 
                         "latitude" = st_coordinates(argos)[1,2],
                         "lag_release_hrs" = round(as.numeric(t),1),
                         "lag_track_end_hrs" = round(as.numeric(t2),1),
                         "release_type" = summary$ReleaseType)

if(!missing(drift_correction)){ 
  
  corrected_release = drift_correct(track_end,argos,N=drift_correction)
  
 
  result[[d]] = cbind(result[[d]],corrected_release)



}

}
}

do.call(rbind,result)

}


# -----------------------------------------------------------------------------
# drift_correct
# -----------------------------------------------------------------------------

# Corrects for tag drift effects between reported release time (or track end for
# floaters) and the first Argos transmissions by back-extrapolating using average
# drift speed and bearing derived from Argos locations over a given time period

drift_correct = function(track_end,argos,N=24){
  
  if(is(argos,'Spatial')) argos = st_as_sf(argos)
  if(!is(argos,'sf')) argos = subset(argos,!is.na(Longitude)) %>% st_as_sf(coords=c('Longitude','Latitude'),crs=4326) 
  
  drift_track = 
    subset(argos,LocationQuality %in% c('1','2','3')) %>%
    arrange(Date) %>%
    mutate(ellapse = as.numeric(difftime(Date,first(Date),unit='hours'))) %>%
    dplyr::filter(row_number() <= which.min(abs(ellapse-N)))
  
  end_pts = 
  dplyr::filter(drift_track,row_number()==1 | row_number() == n())
  
  dist = as.numeric(st_distance(end_pts [1,],end_pts [2,]))
  time = as.numeric(difftime(end_pts$Date[2],end_pts$Date[1],unit='secs'))
  speed_ms = dist/time
  
  bear = as.numeric(lwgeom::st_geod_azimuth(end_pts$geometry) * (180/pi))
  bear = ifelse(bear<0,360+bear,bear)
  backbear = ifelse(bear>=180,bear-180,bear+180)
  
  d0 = as.numeric(difftime(end_pts $Date[1],track_end,unit='secs')) * speed_ms/1000
  
  xy0 = 
  maptools::gcDestination(st_coordinates(end_pts)[1,1],
                          st_coordinates(end_pts)[1,2],
                          backbear,
                          d0) %>% data.frame() %>%
                          setNames(c('drift_correct_long','drift_correct_lat'))

  print(
  ggplot(drift_track) + 
    geom_sf() + 
    annotate('line',x=st_coordinates(end_pts)[,1],
             y=st_coordinates(end_pts)[,2],colour='blue',arrow=arrow(ends='first')) +
    geom_sf(data =  st_as_sf(xy0,coords=c(1,2),crs=4326),colour='red') + 
    annotate('line',x=c(st_coordinates(end_pts)[1,1],xy0[1,1]),
             y=c(st_coordinates(end_pts)[1,2],xy0[1,2]),colour='red',arrow=arrow())
  )
  
  return(xy0)
  
  
  }


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# gpe3_spot_locs
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#' This function extracts Argos and Fastloc locations from SPOT tags and prints them in a
#' format for batch pasting into the metadata of PSAT tags in WC portal. 
#' All locations after an optional release.date (in POSIXct) are excluded

gpe3_spot_locs = function(dir,release.date=NULL){
  
  files = list.files(dir,pattern="\\Locations.csv$",full.names=T)
  if(length(files)==0) stop ('No Locations.csv file in dir')
  files = files[which.max(nchar(files))] #choose file with most characters as this is Fastloc output
  locs  = read.csv(files,stringsAsFactors = F)
  locs  = locs[locs$Type!='User',]
  dt    = do.call(rbind,lapply(strsplit(locs$Date," "),function(x) {x[1] = gsub("\\..*","",x[1]); x}))
  dtpos = as.POSIXct(paste(dt[,1],dt[,2]),format="%H:%M:%S %d-%b-%Y",tz='GMT')
  out = cbind(dt[,2:1],locs[,c('Latitude','Longitude')])
  if(!is.null(release.date)) out = out[dtpos<release.date,]
  out = out[!duplicated(out), ]
  out = apply(out,2,trimws)
  writeLines(apply(out,1,paste,collapse=","))
  
}


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# read_spot
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# reads a table of SPOT tag locations in from a WC portal download folder, or series of folders. 
# Incorporates both Argos and FastGPS positions in a custom table that includes all relevant location error 
# information, sorts by id and date, formats dates and removes duplicate rows.

read_spot = function(dirs){
  
  # argos
  
  argos = lapply(dirs,function(x){files = list.files(x,pattern="\\Locations.csv$",full.names=T)
  if(length(files)!=0) read.csv(files[which.min(nchar(files))],stringsAsFactors = F)})
  argos[sapply(argos, is.null)]  = NULL
  argos = do.call(rbind,argos)[,c(2,4:12)]
  argos = argos[argos$Type=="Argos",]
  names(argos)[c(2,4)] = c('datetime','lc')
  names(argos) = tolower(names(argos))
  argos[,c('satellites','residual')] = NA
  argos = argos[,c(1:3,5:6,4,11:12,7:10)]
  
  # fastgps
  
  fastloc = lapply(dirs,function(x){files = list.files(x,pattern="\\FastGPS.csv$",full.names=T)
  if(length(files)!=0) read.csv(files[which.max(nchar(files))],skip=3,stringsAsFactors = F)})
  
  fastloc[sapply(fastloc, is.null)]  = NULL
  fastloc = do.call(rbind,fastloc)
  fastloc = fastloc[!is.na(fastloc$Latitude) & !is.na(fastloc$Longitude),c(1:3,9,14:15,18)]
  fastloc$datetime = paste(fastloc$Time,fastloc$Day)
  fastloc[,c('Day','Time')] = NULL
  names(fastloc)[1] = 'ptt'
  fastloc$type = "FastGPS"
  fastloc[,c('lc','error.radius','error.semi.major.axis','error.semi.minor.axis','error.ellipse.orientation')] = NA
  names(fastloc) = tolower(names(fastloc))
  fastloc = fastloc[,c(1,6:7,3:4,8,2,5,9:12)]
  
  # combine
  
  locs = rbind(argos,fastloc)
  
  # format dates
  
  dt = do.call(rbind,lapply(strsplit(locs$datetime," "),function(x) {x[1] = gsub("\\..*","",x[1]); x}))
  locs$datetime = as.POSIXct(paste(dt[,1],dt[,2]),format="%H:%M:%S %d-%b-%Y",tz='GMT')
  
  # remove duplicates
  
  locs = locs[!duplicated(locs),]
  
  # order
  
  locs = locs[order(locs$ptt,locs$datetime),]
  
  
}
