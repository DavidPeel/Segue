ExtractASegue <-  function(segueObj, stratum=NULL, along.track.dist.m)
{

  segname<-paste0('seg',along.track.dist.m,'m')

  segueObj_el <- segueObj[[segname]]

  if (!is.null(stratum))
  {
    for (i in 1:length(stratum))
    {
      segueObj_el$segment.data <- segueObj_el$segment.data[segueObj_el$segment.data[, names(stratum)[i]]==stratum[i],]
      segueObj_el$observation.data <- segueObj_el$observation.data[segueObj_el$observation.data[, names(stratum)[i]]==stratum[i],]
    }
  }

  return(segueObj_el)
}


fn.GIS_HeadAlongBearing <-
  function(lon1,lat1, bearing, distance, rEarth=6371.01)
  {
    radianBearing<-bearing/180*pi*-1
    radialDistance<- distance / rEarth
    lat1<-lat1/180*pi
    lon1<-lon1/180*pi

    lat = asin(sin(lat1)*cos(radialDistance)+cos(lat1)*sin(radialDistance)* cos(radianBearing))

    #     if(cos(lat) == 0) {  # Endpoint a pole
    #        lon=lon1;
    #     }
    #     else {
    #        lon = ((lon1-asin(sin(radianBearing)*sin(radialDistance)/cos(lat))+pi) %% (2*pi)) - pi
    #     }


    lon=((lon1-asin(sin(radianBearing)*sin(radialDistance)/cos(lat))+pi) %% (2*pi)) - pi
    lon[cos(lat)==0]<-lon1[cos(lat)==0]

    lon<-lon/pi*180
    lat<-lat/pi*180

    ans<-cbind(lon,lat)

    ans

  }
