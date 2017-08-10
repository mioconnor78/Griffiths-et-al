### some general biodiv data code


# data prep (melt, etc) ---------------------------------------------------
## Make a wide datafile, if needed: melt and recast so the datafile goes from long to wide (species as columns)
data.m <- melt(data, id = c(1,2,3,4,5,52, 53))

## some cleaning of species names
levels(data.m$variable)[levels(data.m$variable)== "Bittium.spp."] <- "Lirobittium.spp."
levels(data.m$variable)[levels(data.m$variable)== "Olivella.sp."] <- "Callianax.sp."
levels(data.m$variable)[levels(data.m$variable)== "Cypricercus."] <- "Cyprideis.beaconensis"
levels(data.m$variable)[levels(data.m$variable)== "Odontosyllis"] <- "Polychaete1"

## replace NAs with 0s, if desired
data[is.na(data)] <- 0


## if you don't want size class information: 
## sum across size classes within plots (samples)

data.p <- ddply(data, .(site, date1, Sample, Time.Code2, variable, dfw,order.dfw,area,salinity, shoot.density, fetch.jc), summarise, sum(value))

data.p$time.ID <- paste(data.p$site, data.p$Time.Code2, sep = '.') #could look at finer time resolution by using Time.Code here
names(data.p) <- c("site", "Date", "Sample", "Time.Code2", "species", "dfw","order","area","salinity","shoot.density","fetch","abundance", "TimeID")
