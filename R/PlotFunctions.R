plot.segue <- function(segueObj,stratum=NULL, along.track.dist.m  , lwd=3, xlab='Longitude', ylab='Latitude', colmap=rainbow, ncols=12, plotattribute=NULL, addsightings=T, ...)
{


  segueObj_el <- ExtractASegue(segueObj, stratum=stratum, along.track.dist.m)

  segdata <- segueObj_el$segment.data # Just to save typing

  plot(range(segdata$lon.start, segdata$lon.end) , range(segdata$lat.start, segdata$lat.end), type='n', xlab=xlab, ylab=ylab, ...)

  if (is.null(plotattribute))
  {
   cols<-sample(colmap(nrow(segdata)),nrow(segdata), replace = F)
  }
  else
  {
    id<-segdata[,plotattribute]
    id[!is.na(id)]<-cut(id[!is.na(id)], breaks=seq(min(id,na.rm=T),max(id,na.rm=T),length.out=ncols), labels=FALSE)
    id[is.na(id)]<-0
    id<-id+1
    cols<-c('gray',colmap(ncols))[id]
  }

  segments(segdata$lon.start, segdata$lat.start, segdata$lon.end, segdata$lat.end, col=cols, lwd=lwd)

  if (addsightings)
  {
    points(segueObj$rawdata$sightings$Longitude, segueObj$rawdata$sightings$Latitude, pch=19)
    tem<-fn.GIS_HeadAlongBearing(segueObj$rawdata$sightings$Longitude, segueObj$rawdata$sightings$Latitude, bearing=segueObj$rawdata$sightings$Angle, distance=segueObj$rawdata$sightings$distance/1000)
    segments(segueObj$rawdata$sightings$Longitude, segueObj$rawdata$sightings$Latitude, tem[,1], tem[,2], col='grey',lwd=2)
    points(tem[,1], tem[,2], pch='x',col='blue')
  }

}
