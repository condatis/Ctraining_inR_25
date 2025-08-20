###################### Functions to run Condatis in R ##############
##Jenny Hodgson, 2024, with important contributions by Tom Travers and Claudia Gutierrez-Arrellano
##Not for open sharing yet, contact jenny.hodgson@liverpool.ac.uk with any queries

##library dependencies - not run to make it easier to 'source' this file
#library(sf)
#library(terra)
#library(sfheaders)
#library(dplyr)
###

makehabitatpoints<-function(habitatraster,inm=TRUE){
  #converts a habitat raster to points in preparation for the Condatis from points function
  habpt <- terra::as.points(habitatraster, values = TRUE, na.rm = TRUE)
  habptcoord <- terra::crds(habpt, df = TRUE)
  habpt <- as.data.frame(habpt)#or would it be better to use st drop geometry?
  habpt <- cbind(habptcoord, habpt)
  
  if(inm){
    #raster units will tend to be in m, and for Condatis we want km distances
    names(habpt)<- c("xm","ym","cover")
    habpt$x <- (habpt$xm)/1000
    habpt$y <- (habpt$ym)/1000
  }else{
    #if not in m, original co-ordinates will be retained - have not checked on implications of this later on
    names(habpt)<- c("x","y","cover")
  }
  habpt<-habpt[habpt$cover>0,]#remove cells with 0 cover
  habpt
}


stargetpt <- function(straster,inm=TRUE){
  #converts source-target raster, coded 1 and 2, to points in preparation for the Condatis from points function
  pts <- terra::as.points(straster, values=TRUE, na.rm=TRUE)
  rascoord <- terra::crds(straster, df=TRUE)
  stpt1 <- as.data.frame(pts)#or would it be better to use st drop geometry?
  stpt1<-cbind(rascoord, stpt1)

  if(inm){
    #raster units will tend to be in m, and for Condatis we want km distances
    names(stpt1)<- c("xm","ym","st")
    stpt1$x <- (stpt1$xm)/1000
    stpt1$y <- (stpt1$ym)/1000
  }else{
    #if not in m, original co-ordinates will be retained - have not checked on implications of this later on
    names(stpt1)<- c("x","y","st")
  }
  
  targetpt<-stpt1[stpt1$st==2,]
  sourcept<-stpt1[stpt1$st==1,]
  st_list<-list(sourcept, targetpt)
  names(st_list)<-c('sourcept','targetpt')
  return(st_list)
}
#################### core condatis function with several output options, but no GIS outputs ##############

condatis_from_points <- function(habpt, sourcept, targetpt, R, 
                        disp, cellside = 1, metric ='ALL',  filename = NA, 
                        threshold = 0.995, minscore=0.2, maxlink = 50000, minlink = 1, maxdisp = Inf){
  ##if you don't need bottlenecks, use metric="flow" to save time
  ##if filename is NA, no output files will be written but the result will be an R list object, which could later be saved
  ##filename can incorporate a path; it should be the first part of your job identifier, additional labels and endings will be added
  ##the settings threshold, minscore, maxlink and minlink are all to do with the number of potential bottlenecks that are saved
  ##be aware that the total number of links in the circuit is number of cells squared - if you save too many you risk an object too big to assign
  ##for England 2024 work we used threshold=0.995 , minscore=0.2 ,maxlink=50000, minlink=1 - these affect the number of bottlenecks that are saved
  ##deliberately led to more bottlenecks than we needed, to allow for some troubleshooting without having to re-run everything
  ##maxdisp is an optional cutoff for dispersal impossibility, in km
  
  # Get the distances between each cell of habitat and every other cell
  len <- dim( habpt)[1]
  dm <- dist( habpt[, c('x','y')])
  
  origin <- sourcept
  target <- targetpt
  
  # Define alpha (mean dispersal) and normalisation so the area under the dispersal kernel integrates to 1
  alpha <- 2/disp
  norm <- R*alpha^2/2/pi*cellside^4
  
  #### Core Condatis Calculations ---------
  
  #Current between cells
  Cfree <- norm*outer( habpt$cover,  habpt$cover, '*')*exp(-alpha*as.matrix(dm))
  diag(Cfree) <- 0
  
  if(maxdisp!=Inf){#NEW optional ability to make connections beyond a certain distance 0
    Cfree[as.matrix(dm) >= maxdisp] <- 0
  }
  
  #Current into a cell and out of a cell
  Cin <- norm* habpt$cover*
    rowSums(exp(-alpha*sqrt(outer( habpt[, 'x'], origin[, 'x'], '-')^2 +
                              outer( habpt[, 'y'], origin[, 'y'], '-')^2)))
  
  Cout <- norm* habpt$cover*
    rowSums(exp(-alpha*sqrt(outer( habpt[, 'x'], target[, 'x'], '-')^2 +
                              outer( habpt[, 'y'], target[, 'y'], '-')^2)))
  M0 <- diag(Cin + Cout + rowSums(Cfree)) - Cfree
  w <- Cin - Cout
  
  v0 <- solve(M0, w, tol = exp(-255)) # This produces the resistance values in the network and is where errors will most likely occur, see 'try' function for error handling
  
  
  I0 <- (v0 + 1) * Cout
  I1 <- (1 - v0) * Cin
  
  ## Conductance value of whole network ----
  cond <- (sum(I0) + sum(I1))/4 #4 to make overall voltage difference 1 rather than 2
  
  # list of overall statistics for saving
  overall<-list(job=filename,ncells=len,conductance=cond)
  
  ## Flow and progress by cell ------
  flo <- apply(Cfree * outer(v0, v0, '-'), 1, function(x) {
    sum(abs(x))
  })
  
  flo <- ((flo)/2 + I0 + I1)/2
  progress<- (v0 + 1)/2
  
  overall$minflo<-min(flo)
  overall$maxflo<-max(flo)  
  overall$minprog<-min(progress)
  overall$maxprog<-max(progress)
  
  # combine progress, flow, standardised flow, and conductance into a data.frame for saving later
  f <- cbind(habpt, flow = flo, progress = progress,  std_flow = flo/max(flo))
  
  
  
  if(metric == 'flow'){
    
    results <- list(overall,f)
    names(results) <- c('overall','flow')  
    
    if(!is.na(filename)){
      write.csv(f, paste0(filename,'_flow.csv'))
      write.csv(as.data.frame(overall), paste0(filename,'_overall.csv'))
      save(results,file=paste0(filename,'_results.Rdata'))
    }
    
    return(results)
    
    stop()
  }#end if metric is flow
  
  ## Power calculations - for bottlenecks------
  
  powr <- Cfree * (outer(v0, v0, '-')/2)^2 #NEW divide by 2 to get same units again
  powlong <- data.frame(
    a = c(matrix(1:len, nrow = len, ncol = len)[upper.tri(powr)]),
    b = c(matrix(1:len, nrow = len, ncol = len, byrow = T)[upper.tri(powr)]),
    powr=c(powr[upper.tri(powr)]))
  
  powlong <- powlong[order(-powlong$powr), ]#sorting the data frame so highest power comes first
  sumpow <- sum(powlong$powr)#total power
  powlong$thresh <- cumsum(powlong$powr)/sumpow
  
  #new (2023) score of bottleneck
  powlong$score<- powlong$powr/sumpow*len
  overall$sumpow<-sumpow
  
  # subset the power scores that account for threshold of the flow with a minimum of minlink and a maximum of maxlink of the highest powers
  upto1 <- 1+sum( powlong$thresh <= threshold)#threshold based on cumulative proportion of circuit's total power
  upto2 <- sum(powlong$score >= minscore)#NEW condition based on score (minimum level can be chosen)
  upto <- max(upto1,upto2)#NEW the larger of the two limits will be chosen
  
  if ( upto < minlink){#additionally, minlink and maxlink will override any scores based on the power
    powlong <- powlong[1:minlink,]
  } else {
    if (upto > maxlink) {
      powlong <- powlong[1:maxlink,]
    } else{powlong<-powlong[1:upto,]
    }
  }
  
  #clean up to save memory
  rm(Cfree)
  rm(powr)
  gc()  
  
  
  ##### #Create dataframes of power scores and location of ends of the bottleneck
  #powpoints is to convert to line geometries, 2 rows per bottleneck
  #power is one row per bottleneck - to continue analysis in R  
  
  powlong$label <- paste(powlong$a, powlong$b, sep = '_')
  powlong$perc<-powlong$powr/sumpow*100 #percentage of total  
  #overall information on no. bottlenecks
  overall$nbottlenecks<- dim(powlong)[1]
  overall$nabovethreshold<-upto1
  overall$naboveminscore<-upto2
  overall$maxperc<-max(powlong$perc)  
  
  powpoints <- cbind( habpt[powlong$a, c('xm', 'ym')], powlong[,c('label','powr','perc', 'score')], type = 'a')
  powpoints <- rbind(powpoints, cbind( habpt[powlong$b, c('xm', 'ym')], powlong[, c('label','powr','perc','score')], type = 'b'))  

  #extract the co-ordinates of bottleneck ends for saving with the power information
  powxy <- habpt[powlong$a, c('xm', 'ym')]
  names(powxy) <- c('xma', 'yma')
  powxy <- cbind(powxy, habpt[powlong$b, c('xm', 'ym')])
  names(powxy)<-c('xma', 'yma','xmb','ymb')
  
  power <- cbind(powxy,powlong)
  
  #new - calculate bottleneck length in case useful
  power$lengthkm <- with(power,sqrt((xma-xmb)^2+(yma-ymb)^2)/1000)
  
  powpoints <- powpoints[order(powpoints$label),] # the ordering is important for the conversion to shapefile
  
  
  ##These outputs will be implemented unless metric= "flow" (in which case stops at line 81)---
  
  results <- list(overall=overall,flow=f,power=power, bottleneck_ends=powpoints)
  
  if(!is.na(filename)){
    write.csv(f, paste0(filename,'_flow.csv'))
    write.csv(power, paste0(filename,'_power.csv'))
    write.csv(powpoints, paste0(filename,'_bottleneck_ends.csv'))
    write.csv(as.data.frame(overall), paste0(filename,'_overall.csv'))
    
    
    save(results,file = paste0(filename,'_results.Rdata'))
  }
  
  return(results)
  
}#end function


########## function to convert bottleneck endpoints to a line sf object #####################

ends_to_lines<-function(endstab,minscore=0,crs=NA){
  #converts the endpoints given in endstab to a line sf object
  #within endstab, x co-ordinates must be named 'xm', y co-ordinates must be named 'ym', points to be connected into
  #the same line must be denoted by 'label', and the start vs the end of the line should be denoted by 'type'
  #any supplied crs will be assigned to the sf object which is output. 
  #if endstab has a column called 'score' only scores above minscore will result in lines
  #

  if("score" %in% names(endstab)){
ends<-ends[ends$score>=minscore,]#changed to minor
}

lineobj <- sfheaders::sf_linestring(
  obj = endstab,
  x = 'xm',
  y = 'ym',
  z = NULL,
  m = NULL,
  linestring_id = 'label',
  keep = TRUE)

st_crs(lineobj)<-crs#new, set correct crs
lineobj <- subset(lineobj, select = -c(type))#remove the type column, used to distinguish beginning from end of line
lineobj$length <- sf::st_length(lineobj)# calculate length
return(lineobj)
}

###################### NEW function for bottleneck zone polygons given score range, based on Claudia's 2023 work ######

zones_by_score_range <- function(lineobj, score_range,directioncol="direction",regioncol="region",keepcols=c()){ 
  
  #lineobj should be a simple features object, it must have a column called 'score', and one called 'length'
  #score_range should have two elements,the lower and upper bounds of scores to be retained
  #if the names given in directioncol and regioncol match names in lineobj, different text in the input will be concatenated to produce
  #output columns called 'directions' and 'regions' respectively. If a name doesn't match, the step will be ignored
  #any columns named in 'keepcols' will make sense if they don't vary within a zone. If they do vary, only the first instance will be kept
  #when duplicates are removed. Use c() if you don't wish to keep any columns beyond the default ones
  
  `%>%` <- dplyr::`%>%`
   mycrs<- st_crs(lineobj)
  
  b <- dplyr::filter(lineobj, score >= score_range[1] & score < score_range[2])
  if(dim(b)[1]>0){
    #create a layer of points at the middle of bottlenecks
    b_point<- sf::st_line_sample(b, sample=0.5)%>% #allows to add columns
      sf::st_sf()
    
    #create a buffer around the points, with the original attributes
    
    b_buffer<- cbind(b_point,st_drop_geometry(b))
    b_buffer$buf_length<-b_buffer$length/2#calculate  bottleneck mid-length
    
    b_buffer<- sf::st_sf(b_buffer)
    
    b_buffer<-sf::st_buffer(b_buffer, dist=st_drop_geometry(b_buffer$buf_length)) #use mid-length as buffer distance
    
    #dissolve overlapping buffers into irregular polygons (units)
    
    b_units<-b_buffer%>%       
      sf::st_union()%>%
      sf::st_cast('POLYGON')%>%
      sf::st_sf() %>%
      dplyr::mutate(
        unit = dplyr::row_number()) #assigns each polygon an ID
    
    #add information about unit, based on constituent buffers (this involves duplicating the units, then later deleting duplicates)
    
    b_score_sum<- sf::st_join(b_units, b_buffer) %>%
      dplyr::group_by(unit)%>%
      dplyr::mutate(sumscore=sum(score),
                    line_count=length(score) )#calculate the sum of score in each unit, and number of buffers that form each unit
  
    #which columns to keep in output
    keepcols<-c(keepcols,"sumscore","line_count")
    
    #record unique directions and regions, if present    
    if(directioncol %in% names(b)){
      b_score_sum<- mutate( b_score_sum,directions=
                              
                              paste(unique(get(directioncol)),collapse="."))
      keepcols<-c(keepcols,"directions")    #we'll keep the concatenated column not the original column
    }
    
    if(regioncol %in% names(b)){
      b_score_sum<- mutate( b_score_sum,regions=
                              
                              paste(unique(get(regioncol)),collapse="."))
      keepcols<-c(keepcols,"regions")
      
    }
    
    #remove all duplicated records, keeps only units, score sum and number of buffers per unit
    b_score_sum<- b_score_sum[!duplicated(b_score_sum$unit),]
    
    
    #remove unnecessary information, keeping only keepcols, as updated above
    b_score_sum<-b_score_sum[,keepcols]
    
    st_crs(b_score_sum)<-mycrs
    return(b_score_sum)
  }#end if
  return(NA)#if there are no bottlenecks within the supplied score range
}