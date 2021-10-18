dsm_a_segueObj <- function(segueObj, stratum=NULL, along.track.dist.m  ,
                           formula, ddf.obj, engine = "gam",
                           convert.units = 1, family = quasipoisson(link = "log"), group = FALSE,
                           control = list(keepData = TRUE), availability = 1, strip.width = NULL,
                           segment.area = NULL, weights = NULL, transect = "line", method = "REML",
                           ...)
{
  library(dsm)
  library(mvbutils)

  # DP could use ... instead of specifying all dsm arguments?


  segueObj_el <- ExtractASegue(segueObj, stratum=stratum, along.track.dist.m)


  modelfit <- dsm(formula=formula, ddf.obj=ddf.obj, engine = engine,
                  convert.units = convert.units, family = family, group = group,
                  control = control , availability = availability , strip.width = strip.width,
                  segment.area = segment.area , weights = weights, transect = transect, method =  method,
                  segment.data = segueObj_el$segment.data, observation.data = segueObj_el$observation.data, ... )


  return(modelfit)
}
