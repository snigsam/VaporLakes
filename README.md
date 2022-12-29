# VaporLakes
Code to track atmospheric vapor lakes that cross east African coastlines, defined by a CWV value (contour). 

Next, for each lake, gather gridded data in "buffer" shapefiles, defined by distances from the lake boundary. Useful package is recommended: rioxarray, example is here https://gis.stackexchange.com/questions/328128/extracting-data-within-geometry-shape/328320#328320. 

I don't think I have to go all the way to OSGEO packages to define a "mask" for the array, as exampled here https://gis.stackexchange.com/questions/354782/masking-netcdf-time-series-data-from-shapefile-using-python/354798#354798


