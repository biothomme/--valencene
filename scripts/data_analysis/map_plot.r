devtools::install_github("joeroe/rpaleoclim")
install.packages('rgdal')

library("rgdal")
library("rpaleoclim")
library("raster")

# Bio_1=Annual Mean Temperature [°C*10]
# Bio_2=Mean Diurnal Range [°C]
# Bio_3=Isothermality [Bio_2/Bio_7]
# Bio_4=Temperature Seasonality [standard deviation*100]
# Bio_5=Max Temperature of Warmest Month [°C*10]
# Bio_6=Min Temperature of Coldest Month [°C*10]
# Bio_7=Temperature Annual Range [°C*10]
# Bio_8=Mean Temperature of Wettest Quarter [°C*10]
# Bio_9=Mean Temperature of Driest Quarter [°C*10]
# Bio_10=Mean Temperature of Warmest Quarter [°C*10]
# Bio_11=Mean Temperature of Coldest Quarter [°C*10]
# Bio_12=Annual Precipitation [mm/year]
# Bio_13=Precipitation of Wettest Month [mm/month]
# Bio_14=Precipitation of Driest Month [mm/month]
# Bio_15=Precipitation Seasonality [coefficient of variation]
# Bio_16=Precipitation of Wettest Quarter [mm/quarter]
# Bio_17=Precipitation of Driest Quarter [mm/quarter]
# Bio_18=Precipitation of Warmest Quarter [mm/quarter]
# Bio_19=Precipitation of Coldest Quarter [mm/quarter]

lh_10m <- paleoclim("lh", "10m", region = extent(c(-15, 50, 30, 50)))
annual_temp = lh_10m[[c(1)]]
annual_temp =  setMinMax(annual_temp)

cellStats(annual_temp, min)
cellStats(annual_temp, max)

image(annual_temp)
col <- terrain.colors(7)[c(2:4)]
brk <- c(-60,200, 265)

plot(annual_temp, col=col, breaks= brk)
