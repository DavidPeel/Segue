# Code to do segmentation of data ready for dsm



MakeSegue <- function(sightings, effort, trunc.dist=3000, truncate_sdata=T, along.track.dist.m=c(2000,5000,10000),
                        strata.cols =  c(), enviro.covars = c(), lon.col = "Lon", lat.col = "Lat", effort.datetime.col = "date.time",   #strata.cols = c("Area","Migration"),
                        sighting.no.col = "sighting.no.unique", sight.datetime.col = "Time", transect.sub.name = "TransectID",size.name ="Count", xyprojection=NULL)
{
  #
  #Input/arguments
  # sightings = name of the original sighting data frame
  # effort = name of the original effort data frame
  # trunc.dist = what is the truncation distance for the fitted distance model? (in units of 'distance' column given in sightings data frame)
  # truncate_sdata = T/F if you want to truncate the sighting data
  # along.track.dist.m is the along track distance you'd like to segment into (can be a vector of different distances to do)
  # strata.cols =
  # enviro.covars = a LIST of column names which contains the environmental/sightability covars that remained in the detection function AFTER det func model selection has been undertaken.
  # if you'd like strata/region, etc, information to come through into the final segment dataframe, name them here as they appear as column names
  # lon.col = which column (as a text string) has the longitude data
  # lat.col = which column (as a text string) has the latitude data
  # datetime.col = which column (as a text string) has the date time object data
  #
  #
  # transect.sub.name = which column (as a text string), which has the continous chunks of effort. This could be the entire transect, if there were no 'off effort' bits
  #                     within a single transect. Otherwise, you'll have to have gone through the effort to work out unique chunks of effort time that were either on or off effort.
  # size.name = which column (as a text string) has the group size information for each sighting?
  # xyprojection = The R projection string if you want to work with projection



  # Preliminaries -----------------------------------------------------------

  library( geosphere)
  library( zoo)
  library( stringr)
  library(tidyr)


  # Process DATA -----------------------------------------------------------


  if (truncate_sdata)
  {
    sightings <- sightings[sightings$distance<=trunc.dist, ]
  }


  #Modify 'sightings' dataframe
  sightings$sighting.no.unique<- 1:dim( sightings)[1]

  sightings$Transect.alt <- paste( sightings$Transect, sightings$Migration, sep="_") #create a new transect name that will appear in both sightings and effort data, and one that separates out north and south migration


  #Modify 'effort' dataframe

  effort$effort.status<- rep(1, dim(effort)[1]) #below code requires an 'effort.status' column; 1 == on effort, 0 == off effort. Here there is no off effort

  effort.2rows<- do.call( rbind, lapply( 1:dim(effort)[1], conv.to.2.rows.f, effort.data = effort))

  effort.2rows$Transect.old<- effort.2rows$Transect #relabel the Transect column
  effort.2rows$Transect<- effort.2rows$TransectID  #relabel the TransectID to be the new Transect column

  effort.2rows$Transect.alt<- paste( effort.2rows$Transect.old, effort.2rows$Migration, sep="_") #create a new transect name that will appear in both sightings and effort data, and one that separates out north and south migration


  #ans<-list(segment.data=list(), observation.data=list(), rawdata=list())
  ans<-list()

  segment.summaries.all<-list()
  unique.segment.no<-list()
  summarise.sightings.at <- list()
  summaries <- list()
  obsdata <-list()

  for (idist in along.track.dist.m)
  {
    iname <- paste0('seg',idist,'m')

    ans[[iname]] <- list()

    # TRANSECT SEGMENTATION -----------------------------------------------------------
    #Run the 'transect.segmentation.f' function to segment the effort data into predefined lengths
    segment.summaries.all <- transect.segmentation.f( effort.dataset = effort.2rows, transect.sub.name = transect.sub.name , lon.col = lon.col, lat.col = lat.col,
                                                              datetime.col = effort.datetime.col, along.track.dist.m = idist, enviro.covars = enviro.covars, strata.cols = strata.cols)

    # MATCHING SIGHTINGS AND EFFORT -----------------------------------------------------------
    #match sightings at  segment level
    unique.segment.no <- unlist(lapply( sightings$sighting.no.unique, segment.no.per.sighting.f, segment.summary = segment.summaries.all,
                                                sighting.data = sightings, sighting.no.col = sighting.no.col, date.time.col = sight.datetime.col))
    datacolm <- paste0('unique.segment.no.',idist,'m')
    sightings[ ,datacolm]<- unique.segment.no

    # SUMMARISE SEGMENTS -----------------------------------------------------------
    summarise.sightings.at <- summarise.sightings.at.segment.f( segment.data = segment.summaries.all, sighting.summary.data = sightings,
                                                                        sighting.summary.data.column = datacolm, trunc.dist = trunc.dist, size.name = size.name)

    ans[[iname]]$segment.data <- cbind( segment.summaries.all, summarise.sightings.at)


    # Hack to rename a couple of columns so it matches what DSM expects - DP
    colnames(ans[[iname]]$segment.data)[match("unique.segment.no",colnames(ans[[iname]]$segment.data))] <- "Sample.Label"
    colnames(ans[[iname]]$segment.data)[match("segment.length.m",colnames(ans[[iname]]$segment.data))] <- "Effort"

    if (!is.null(xyprojection))
    {
      tem<- fn.projection(ans[[iname]]$segment.data$lon.mid, ans[[iname]]$segment.data$lat.mid,ToPrj = xyprojection)
      ans[[iname]]$segment.data$x <- tem[[1]]
      ans[[iname]]$segment.data$y <- tem[[2]]
    }
    else
    {
      ans[[iname]]$segment.data$x <-ans[[iname]]$segment.data$lon.mid
      ans[[iname]]$segment.data$y <- ans[[iname]]$segment.data$lat.mid
    }

    # Link seg data to sightdata to make obs data
    id <- match(sightings[ , datacolm], ans[[iname]]$segment.data$Sample.Label)


    ans[[iname]]$observation.data<- data.frame(object=sightings$sighting.no.unique, Sample.Label=sightings[ , datacolm],
                                   size=sightings$Count, distance=sightings$distance, Effort=ans[[iname]]$segment.data$Effort[id], Area=sightings$Area, Migration=sightings$Migration, stringsAsFactors = F)

  }

  # Not sure whether we need this data ever? Useful for plotting I guess
  ans$rawdata=list(sightings=sightings, effort= effort.2rows)

  return(ans)
}


segment.no.per.sighting.f<- function( sighting.no, segment.summary, sighting.data, sighting.no.col, date.time.col)
{
  #print(sighting.no)

  ##Function does:
  # Matches sightings to the unique segments, as defined in the function 'transect.segmentation.f'

  ##Arguments

  # sighting.no = unique sighting number for each individual sighting
  # segment.summary = what is the name of the segment summary object?
  # sighting.data = name of the original sighting data frame
  # sighting.no.col = which column (as a text string) contains the unique sighting number for each individual sighting?
  # date.time.col = which column (as a text string) has the date time object data for each individual sighting?

  ## Returns
  # a vector (of length of the total number of sightings) containing the unique segment number for each sighting


  sighting.data.i<- sighting.data[ which( sighting.data[, match(sighting.no.col, colnames(sighting.data))] == sighting.no),]
  sighting.data.i.date.time<- eval(parse(text= paste("sighting.data.i$",date.time.col,sep="")))


  segment.summary.transect<- segment.summary[ which( segment.summary$Transect.alt ==  sighting.data.i$Transect.alt),]

  #sighting.transect<- sighting.data.i[, match(transect.col, colnames(sighting.data.i))]

  #segment.summary.transect<- segment.summary[ which( segment.summary[, match("global.transect.no", colnames(segment.summary))] == sighting.transect),]


  segment.of.sighting<- findInterval(sighting.data.i.date.time, segment.summary.transect$datetime.start)

  segment.summary.transect.i<- segment.summary.transect[ segment.of.sighting,]

  sighting.unique.segment.no<- segment.summary.transect.i$unique.segment.no

  sighting.unique.segment.no
}




summarise.sightings.at.segment.f<- function(segment.data, sighting.summary.data, sighting.summary.data.column, trunc.dist, size.name)
{
  ##Function does:
  # Calculates the numbers of sightings (and total individuals) per segment of a given length

  ##Arguments
  # segment.data = what is the name of the segment summary object?
  # sighting.summary.data = original sighting data frame
  # sighting.summary.data.column = which column (as a text string) in the sighting data contains the unique segment number for a given segment length
  # trunc.dist = what is the truncation distance for the fitted distance model? (in units of 'distance' column given in sightings data frame)
  # size.name = which column (as a text string) has the group size information for each sighting?

  ## Returns
  # a data frame, with each row give the number of sightings per segment (sightings.n), and the total number of animals per segment (total.n)


  sighting.segment.nos<- eval(parse(text= paste("sighting.summary.data$", sighting.summary.data.column, sep="")))



  sightings.per.segment.f<- function(i)
  {
    sighting.summary.data.seg<- sighting.summary.data[ which( sighting.segment.nos == i),]
    sighting.summary.data.seg<- sighting.summary.data.seg[ which(sighting.summary.data.seg$distance<= trunc.dist),] #make sure all sightings are within truncation distance
    size.col<- eval(parse(text= paste("sighting.summary.data.seg$", size.name, sep="")))

    sightings.n<- dim( sighting.summary.data.seg)[1]
    if( sightings.n == 0) { total.n<- 0} else {total.n<- sum(size.col, na.rm=T) }

    #sightings.calves<- sighting.segment.nos[ which( sighting.segment.nos$calves>0),]
    #sightings.calves.n<- dim( sightings.calves)[1]
    #if( sightings.calves.n == 0) { sightings.calves.total<- 0} else {sightings.calves.total<- sum( sightings.calves$size, na.rm=T)}

    #sightings.no.calves<- sighting.segment.nos[ which( sighting.segment.nos$calves==0),]
    #sightings.no.calves.n<- dim( sightings.no.calves)[1]
    #if( sightings.no.calves.n == 0) { sightings.no.calves.total<- 0} else {sightings.no.calves.total<- sum( sightings.no.calves$size, na.rm=T)}

    #ans.df<- data.frame( i, sightings.n, total.n, sightings.calves.n, sightings.calves.total, sightings.no.calves.n, sightings.no.calves.total)
    ans.df<- data.frame( i, sightings.n, total.n)
    return(ans.df)
  }
  sightings.per.segment<- do.call( rbind, lapply( segment.data$unique.segment.no, sightings.per.segment.f))

  sight.seg.ans.df<- sightings.per.segment
  return( sight.seg.ans.df)
}


transect.segmentation.f<- function( effort.dataset, transect.sub.name, lon.col, lat.col , datetime.col, along.track.dist.m, enviro.covars = NULL, strata.cols )
{
  ##Function does:
  #create 'continuous' segments based on either the pre-defined 'Transect.sub' column and the nominated enviro-covariates, if nominated. If there are no enviro-covariates nominated, then it just uses Transect.sub,

  ##Arguments

  # effort.dataset = effort dataframe
  # transect.sub.name = which column (as a text string), which has the continous chunks of effort. This could be the entire transect, if there were no 'off effort' bits within a single transect. Otherwise, you'll have to have gone through the effort to work out unique chunks of effort time that were either on or off effort.
  # lon.col = which column (as a text string) has the longitude data
  # lat.col = which column (as a text string) has the latitude data
  # datetime.col = which column (as a text string) has the date time object data
  # along.track.dist.m is the along track distance you'd like to segment into
  # enviro.covars = a LIST of column names which contains the environmental/sightability covars that remained in the detection function AFTER det func model selection has been undertaken.
  # if you'd like strata/region, etc, information to come through into the final segment dataframe, name them here as they appear as column names

  ## Returns
  # a dataframe containing a unique row for every new segment.


  if( is.null( enviro.covars)) { transect.sub.col<- eval(parse(text= paste("effort.dataset$", transect.sub.name, sep="")))

  unique.transect.nos<- transect.sub.col
  unique.transect.nos<- str_pad( unique.transect.nos, max(nchar(unique.transect.nos))+1, side = c("left"), pad = "0")


  } else {


    enviro.covars.df<- effort.dataset[,c( transect.sub.name, enviro.covars)]

    #create vector to describe unique along-track combinations of the enviro/sightability covars that remain in the detection function.
    enviro.covars.concat<- unite(data = enviro.covars.df, col = "enviro.covars.concat", colnames( enviro.covars.df), sep = "_")

    rle.lengths<- rle( enviro.covars.concat$enviro.covars.concat)$lengths

    unique.enviros<- rep( 1:length(rle.lengths), rle.lengths)

    unique.transect.nos<- str_pad( unique.enviros, max(nchar(unique.enviros))+1, side = c("left"), pad = "0")}

  effort.dataset$unique.enviros<- unique.transect.nos




  cycle.on.transect.no<- function( unique.enviros.i)
  {
    #print(unique.enviros.i)
    trans.sub.index<- which( effort.dataset$unique.enviros == unique.enviros.i)
    effort.dataset.trans.sub<- effort.dataset[ trans.sub.index,]
    #global.transect.no<- eval(parse(text= paste("effort.dataset.trans.sub$", transect.name, sep="")))[1]

    global.transect.no<- effort.dataset.trans.sub$Transect[1] #get the 'global' transect name
    Transect.alt<- effort.dataset.trans.sub$Transect.alt[1]

    trans.global.index<- which( effort.dataset$Transect == global.transect.no)
    #trans.global.index<- eval( parse( text = paste( "which( effort.dataset$", transect.name, '==global.transect.no', sep="")

    effort.dataset.trans<- effort.dataset[ trans.global.index,]

    transect.sub.index<- which( effort.dataset.trans$unique.enviros == unique.enviros.i)

    global.transect.sub.no<- eval(parse(text= paste("effort.dataset.trans.sub$", transect.sub.name, sep="")))[1]

    if( min(trans.sub.index) == max(trans.global.index)) {NULL} else {

      if( max( trans.sub.index) ==  max(trans.global.index)) {  effort.trans.sub.final<- effort.dataset.trans.sub} else { effort.trans.sub.final<- effort.dataset[ c(trans.sub.index, (max(trans.sub.index)+1)),]}

      effort.dataset.trans<- effort.trans.sub.final



      lon.column<- eval(parse(text= paste("effort.dataset.trans$", lon.col, sep="")))
      lat.column<- eval(parse(text= paste("effort.dataset.trans$", lat.col, sep="")))
      datetime.column<- eval(parse(text= paste("effort.dataset.trans$", datetime.col, sep="")))

      #pull out the existing coordinates (lats/longs and time)
      p<- data.frame(lon.column, lat.column, datetime.column) #all survey coordinates and date/time stamps
      p.coords<- as.matrix( p[,c(1,2)]) #survey waypoints
      p.coords.n<- dim( p.coords)[1]    #number of survey waypoints

      #pull out the pertinent environmental/sightability covariates





      #alongtrack distances for survey waypoints
      along.track.dist.f<- function(i)
      {
        d.i<- distVincentyEllipsoid( p.coords[i,], p.coords[i+1,])
        d.i
      }
      d<- c(0, unlist(lapply( 1:(p.coords.n-1), along.track.dist.f))) # along track distances between each of the successive existing points. In metres

      d.cumsum<- cumsum(d) #cumsum of distances between survey waypoints
      d.total<- sum(d, na.rm=T) #total track length between start and end points. In metres.

      survey.wp.inter.point.dist.m<- c(d.cumsum[1], d.cumsum[2:length(d.cumsum)]- d.cumsum[1:length(d.cumsum)-1]) #estimate the distances between the survey waypoints




      #create cutpoint distances based on total along track distances
      d.cumsum.cutpoints<- c( seq(0, d.total, along.track.dist.m), d.total) # distance is in metres

      #if( all.equal( d.cumsum.cutpoints, survey.wp.inter.point.dist.m) == T) { d.cumsum.cutpoints<- c()} else {d.cumsum.cutpoints<- d.cumsum.cutpoints[2:(length(d.cumsum.cutpoints)-1)]}
      # delete the first and last points of this sequence as these are existing survey waypoints
      if( max(d.cumsum) <= along.track.dist.m) { d.cumsum.cutpoints<- c()} else {d.cumsum.cutpoints<- d.cumsum.cutpoints[2:(length(d.cumsum.cutpoints)-1)]}


      if( length(d.cumsum.cutpoints) == 0) {d.cumsum.cutpoints<- vector()} else {NULL}
      no.points<- length( d.cumsum.cutpoints) #number of individual cutpoints in the segment transect
      cutpoints.index<- rep( 2, no.points)



      dist.df.1<- data.frame( d.cumsum.cutpoints, cutpoints.index)
      lon.cp<- rep(NA, dim(dist.df.1)[1])
      lat.cp<- rep(NA, dim(dist.df.1)[1])
      datetime.cp<- as.POSIXlt( rep(NA, dim(dist.df.1)[1]))
      inter.point.dist.m.NA<- rep(NA, dim(dist.df.1)[1])
      effort.status.NA<- rep(NA, dim(dist.df.1)[1])

      if( length(enviro.covars)==0) {  dist.df.1<- cbind( dist.df.1, lon.cp, lat.cp,datetime.cp, inter.point.dist.m.NA,effort.status.NA)
      colnames( dist.df.1)<- c("d.cumsum", "survey.wp.index", "lon.column","lat.column","datetime.column", "survey.wp.inter.point.dist.m", "effort.status")

      } else{

        enviro.covars.empty.df<- data.frame(matrix(vector(), dim(dist.df.1)[1], length(enviro.covars)))
        colnames( enviro.covars.empty.df)<- paste( enviro.covars, ".NA", sep="")

        dist.df.1<- cbind( dist.df.1, lon.cp, lat.cp,datetime.cp, inter.point.dist.m.NA, enviro.covars.empty.df,effort.status.NA)
        colnames( dist.df.1)<- c("d.cumsum", "survey.wp.index", "lon.column","lat.column","datetime.column", "survey.wp.inter.point.dist.m", enviro.covars, "effort.status")

      }
      #create cutpoint middle coordinates (i.e., a coordinate in the middle of the created segments)

      d.cumsum.cutpoints.mids<- d.cumsum.cutpoints-(along.track.dist.m/2)

      d.cumsum.cutpoints.mids<- c( d.cumsum.cutpoints.mids, sum( c( seq(0, d.total, along.track.dist.m), d.total)[ c( length(c( seq(0, d.total, along.track.dist.m), d.total))-1, length(c( seq(0, d.total, along.track.dist.m), d.total)))], na.rm=T)/2)



      no.points.mid<- length( d.cumsum.cutpoints.mids) #number of individual cutpoints in the segment transect
      cutpoints.mid.index<- rep(3, no.points.mid)
      dist.df.1.mid<- data.frame( d.cumsum.cutpoints.mids, cutpoints.mid.index)
      lon.cp.mid<- rep(NA, dim(dist.df.1.mid)[1])
      lat.cp.mid<- rep(NA, dim(dist.df.1.mid)[1])
      datetime.cp.mid<- as.POSIXlt( rep(NA, dim(dist.df.1.mid)[1]))
      inter.point.dist.m.NA<- rep(NA, dim(dist.df.1.mid)[1])
      effort.status.mid.NA<- rep(NA, dim(dist.df.1.mid)[1])

      #bss.mid.NA<- rep(NA, dim(dist.df.1.mid)[1])
      #cloud.cover.mid.NA<- rep(NA, dim(dist.df.1.mid)[1])
      #turbidity.mid.NA<- rep(NA, dim(dist.df.1.mid)[1])
      if( length(enviro.covars)==0) { dist.df.1.mid<- cbind( dist.df.1.mid, lon.cp.mid, lat.cp.mid,datetime.cp.mid, inter.point.dist.m.NA,effort.status.mid.NA)
      colnames( dist.df.1.mid)<- c("d.cumsum", "survey.wp.index", "lon.column","lat.column","datetime.column", "survey.wp.inter.point.dist.m", "effort.status")

      } else{

        enviro.covars.mid.empty.df<- data.frame(matrix(vector(), dim(dist.df.1.mid)[1], length(enviro.covars)))
        colnames( enviro.covars.mid.empty.df)<- paste( enviro.covars, ".mid", sep="")


        dist.df.1.mid<- cbind( dist.df.1.mid, lon.cp.mid, lat.cp.mid,datetime.cp.mid, inter.point.dist.m.NA, enviro.covars.mid.empty.df,effort.status.mid.NA)
        colnames( dist.df.1.mid)<- c("d.cumsum", "survey.wp.index", "lon.column","lat.column","datetime.column", "survey.wp.inter.point.dist.m", enviro.covars, "effort.status")

      }


      effort.status<- effort.dataset.trans$effort.status
      survey.wp.index<- rep(1, p.coords.n)

      if( length(enviro.covars)==0) {  dist.df.2<- data.frame( d.cumsum, survey.wp.index, lon.column, lat.column, datetime.column, survey.wp.inter.point.dist.m,effort.status )} else{

        enviro.covars.values.df<- as.data.frame( effort.dataset.trans[, match(enviro.covars, colnames(effort.dataset.trans))])
        colnames( enviro.covars.values.df)<- enviro.covars

        dist.df.2<- data.frame( d.cumsum, survey.wp.index, lon.column, lat.column, datetime.column, survey.wp.inter.point.dist.m, enviro.covars.values.df,effort.status )
      }


      if( dim(dist.df.1)[1] == 0) { dist.df.3<- rbind( dist.df.1.mid, dist.df.2)} else {


        dist.df.3<- rbind( dist.df.1, dist.df.1.mid, dist.df.2)
      }

      dist.df.3<- dist.df.3[ order(dist.df.3$d.cumsum),]

      dist.df.3$inter.point.dist.m<- c(dist.df.3$d.cumsum[1], dist.df.3$d.cumsum[2:dim(dist.df.3)[1]]- dist.df.3$d.cumsum[1:(dim(dist.df.3)[1]-1)]) #distances between each ordered survey wp and cutpoints, in metres

      cutpoint.i<- which( dist.df.3$survey.wp.index == 2)
      survey.wp.i<-  which( dist.df.3$survey.wp.index == 1)
      cutpoint.mid.i<- which( dist.df.3$survey.wp.index == 3)
      cutpoints.all<-  sort( c(which( dist.df.3$survey.wp.index == 2) , which( dist.df.3$survey.wp.index == 3)))

      dist.df.4<- dist.df.3  #create dataframe to receive results of for loop below

      dist.df.4$Transect<- effort.dataset.trans$Transect[1]

      #eval(parse(text= paste("effort.dataset.trans$", lon.col, sep="")))

      dist.df.4$Transect.sub<- eval(parse(text= paste("effort.dataset.trans$", transect.sub.name, "[1]", sep="")))


      for(i in cutpoints.all)
      {

        survey.wp.ii<- which( is.na( dist.df.4$lon.column)==F)

        interval.wp<- findInterval( i, survey.wp.ii) #which survey waypoint interval does cutpoint i fall into?

        survey.wp.prior<- dist.df.4[survey.wp.ii[interval.wp],]
        survey.wp.after<- dist.df.4[survey.wp.ii[interval.wp+1],]

        d.survey.wpts<- distVincentyEllipsoid( c(survey.wp.prior$lon.column, survey.wp.prior$lat.column), c(survey.wp.after$lon.column, survey.wp.after$lat.column))

        b.i<- bearing( c(survey.wp.prior$lon.column, survey.wp.prior$lat.column), c(survey.wp.after$lon.column, survey.wp.after$lat.column))
        coord.i<- destPoint(c(survey.wp.prior$lon.column, survey.wp.prior$lat.column), b.i, dist.df.3$inter.point.dist.m[i])

        transect.time.secs.i<- difftime(survey.wp.after$datetime.column, survey.wp.prior$datetime.column, units=c("secs")) #total transect time in seconds

        speed.mean.i<- d.survey.wpts/as.numeric(transect.time.secs.i) # speed in metres per second; mean speed between the two 'transect' points (i.e., start and finish)
        speed.mean.i.inv<- 1/speed.mean.i #mean seconds per metre

        secs.to.point.i<- dist.df.3$inter.point.dist.m[i]* speed.mean.i.inv


        datetime.alongtrack.i<- survey.wp.prior$datetime+secs.to.point.i #estimate the date/time for each point along track (based on average speed). As a strptime object

        dist.df.4$lon.column[i]<- coord.i[1,1]
        dist.df.4$lat.column[i]<- coord.i[1,2]
        dist.df.4$datetime.column[i]<- datetime.alongtrack.i
        dist.df.4
      }

      #create segment membership vector
      cutpoint.i<- c(1, (which( dist.df.4$survey.wp.index == 2)+1))

      membership.i<- rep(NA, dim(dist.df.4)[1])
      membership.i<- replace( membership.i, cutpoint.i, 1:length( cutpoint.i))

      membership.i<- na.locf(membership.i)

      dist.df.4$membership.i<- membership.i  # an index indicating which n km segment each effort row of a given transect

      ### create membership string



      #pad out transect number to be a specific length in characters
      transect.no.padded<- eval(parse(text= paste("str_pad( dist.df.4$Transect.sub, (max( nchar(effort.dataset$", transect.sub.name, "))+1), side = c('left'), pad = '0')", sep="")))

      #pad out membership number to be a specific length in characters.
      membership.no.padded<- str_pad(as.character( dist.df.4$membership.i), max( nchar(as.character( dist.df.4$membership.i)))+1, side = c("left"), pad = "0")


      # add new unique segment number (unique to a given survey, which already has unique transect names)
      dist.df.4$unique.segment.no<- paste( transect.no.padded, "_", membership.no.padded, sep="")

      dist.df.4$effort.status<-  na.locf( dist.df.4$effort.status)

      if( length(enviro.covars)==0) {NULL} else {

        enviro.covars.values.df<- na.locf( dist.df.4[, match(enviro.covars, colnames(dist.df.4))])
        dist.df.4[,match(enviro.covars, colnames(dist.df.4))]<- enviro.covars.values.df
      }
      #dist.df.4$bss<- na.locf( dist.df.4$bss)
      #dist.df.4$cloud.cover<- na.locf( dist.df.4$cloud.cover)
      #dist.df.4$turbidity<- na.locf( dist.df.4$turbidity)


      #summarise at the segment level

      segment.summary.f<- function( unique.segment.no)
      {
        #print( unique.segment.no)
        segment.index<- which( dist.df.4$unique.segment.no == unique.segment.no)
        segment.data.i<- dist.df.4[ segment.index,]

        if( segment.data.i$membership.i[1]==1) {

          lon.start<- segment.data.i$lon.column[ 1]
          lat.start<- segment.data.i$lat.column[ 1]
          datetime.start<- segment.data.i$datetime.column[ 1]} else {  lon.start<- dist.df.4$lon.column[ min(segment.index)-1]
          lat.start<- dist.df.4$lat.column[ min(segment.index)-1]
          datetime.start<- dist.df.4$datetime.column[ min(segment.index)-1]}



        lon.mid<- segment.data.i$lon.column[ which( segment.data.i$survey.wp.index == 3)]
        lat.mid<- segment.data.i$lat.column[ which( segment.data.i$survey.wp.index == 3)]
        datetime.mid<-  segment.data.i$datetime.column[ which( segment.data.i$survey.wp.index == 3)]

        lon.end<- segment.data.i$lon.column[ dim(segment.data.i)[1]]
        lat.end<- segment.data.i$lat.column[ dim(segment.data.i)[1]]
        datetime.end<- segment.data.i$datetime.column[ dim(segment.data.i)[1]]

        #esa.km2.i<- segment.data.i$esw*(segment.data.i$inter.point.dist.m/1000) *2 ## multiplied by 2 because effort on both sides of aircraft)

        segment.length.m<- sum(segment.data.i$inter.point.dist.m, na.rm=T)
        #esa.km2<- sum(esa.km2.i, na.rm=T)

        effort.status<-  segment.data.i$effort.status[1]

        if( length(enviro.covars)==0&&length( strata.cols)==0 ) {
          ans.df<-  data.frame( global.transect.no, global.transect.sub.no, unique.segment.no,Transect.alt,  lon.start, lat.start, datetime.start, lon.mid, lat.mid, datetime.mid, lon.end, lat.end, datetime.end, segment.length.m,effort.status)} else if(length(enviro.covars)==0&&length( strata.cols)>0) {

            strata.cols.df<- as.data.frame( effort.dataset.trans[1, match(strata.cols, colnames(effort.dataset.trans))])
            colnames( strata.cols.df)<- strata.cols

            ans.df<-  data.frame( strata.cols.df, global.transect.no, global.transect.sub.no, unique.segment.no, Transect.alt, lon.start, lat.start, datetime.start, lon.mid, lat.mid, datetime.mid, lon.end, lat.end, datetime.end, segment.length.m,effort.status)} else if ( length(enviro.covars)>0&&length( strata.cols)==0) {

              enviro.covars.mean.df<- as.data.frame( segment.data.i[1,match(enviro.covars, colnames(segment.data.i))])
              # DP Hack to sync with DSM
              #colnames( enviro.covars.mean.df)<- paste( enviro.covars, ".mean", sep="")
              colnames( enviro.covars.mean.df)<- enviro.covars
              ans.df<-  data.frame(global.transect.no, global.transect.sub.no, unique.segment.no, Transect.alt, lon.start, lat.start, datetime.start, lon.mid, lat.mid, datetime.mid, lon.end, lat.end, datetime.end, segment.length.m, enviro.covars.mean.df, effort.status)} else{
                enviro.covars.mean.df<- as.data.frame( segment.data.i[1,match(enviro.covars, colnames(segment.data.i))])
                # DP Hack to sync DSM
                # colnames( enviro.covars.mean.df)<- paste( enviro.covars, ".mean", sep="")
                colnames( enviro.covars.mean.df)<-  enviro.covars
                strata.cols.df<- as.data.frame( effort.dataset.trans[1, match(strata.cols, colnames(effort.dataset.trans))])
                colnames( strata.cols.df)<- strata.cols



                ans.df<- data.frame( strata.cols.df, global.transect.no, global.transect.sub.no, unique.segment.no, Transect.alt, lon.start, lat.start, datetime.start, lon.mid, lat.mid, datetime.mid, lon.end, lat.end, datetime.end, segment.length.m, enviro.covars.mean.df,effort.status)}

        return(ans.df)

      }
      segment.summary<- do.call( rbind, lapply( unique(dist.df.4$unique.segment.no), segment.summary.f))


      return(segment.summary)}

  }
  segment.summary.2<- do.call( rbind, lapply( unique(effort.dataset$unique.enviros), cycle.on.transect.no))

  return( segment.summary.2)

}


# a function to takes the effort start/stop columns in the effort data and pop these in separate rows (as the following code requires)
conv.to.2.rows.f<- function( i, effort.data)
{
  effort.data.i<- effort.data[ i,]

  date.time<- c( effort.data.i$TimeOn, effort.data.i$TimeOut)
  Lon<- c( effort.data.i$Lon0, effort.data.i$Lon1)
  Lat<- c( effort.data.i$Lat0, effort.data.i$Lat1)

  effort.data.i<- cbind( effort.data.i[c(1,1),], date.time, Lon, Lat)
  return( effort.data.i)

}


fn.projection <-  function (x, y, FromPrj="+proj=longlat +datum=WGS84", ToPrj="+proj=utm +zone=53 ellps=WGS84", outformat='list')
{
  library(rgdal)

  xy <- data.frame(ID = 1:length(x), X = x, Y = y)
  coordinates(xy) <- c("X", "Y")
  proj4string(xy) <- CRS(FromPrj)   # CRS("+proj=longlat +datum=WGS84")  ## for example
  res <- spTransform(xy, CRS(ToPrj))  # paste("+proj=utm +zone=",zone," ellps=WGS84",sep='')

  if (outformat=='list')
  {
    res <- as.list(as.data.frame(res)[,-1])
  }
  else if (outformat=='data.frame')
  {
    res<- as.data.frame(res)[,-1]
  }
  else if (outformat=='matrix')
  {
    res<- as.matrix(as.data.frame(res)[,-1])
  }
  else if (outformat=="sp")
  {
    # Do nothing
  }
  else
  {
    stop('Unknown output format')
  }


  return(res)
}
