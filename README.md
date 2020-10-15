# CR18
Codes associated with Cattiaux &amp; Ribes (2018), Defining single extreme weather events from a climate perspective, BAMS.

Details:

- rbase.R contains basic functions, e.g. the treatment of netcdf files (function myno() imports a netcdf in a list, etc.) or the treatment of time/dates, etc.
- functions.R contains the functions dedicated to the event definition issue, e.g. compute.stat.1d() that computes p1, p0, far, etc. at a given location for various time windows and compute.stat.2d() that generalizes to various spatial domains.
- figures.R contains the functions used for plotting the results.
- lon-lat_analysis.R contains the generation of spatial domains of various sizes for a lon-lat analysis.
- country_analysis.R contains the generation of spatial domains of various sizes for a country analysis.

- main.R is the executable script that gets the data, performs the analysis and plots the results.
