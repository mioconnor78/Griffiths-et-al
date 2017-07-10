# Zostera by Smithora
# Prepared for Mary, July 10 2017


setwd("/Users/aolson/Desktop/Masters_Data") # change working directory
data<-read.csv("Macrophyte_biomass_zeros_20170709.csv", stringsAsFactors=F) 

smith<- subset(data,data$macrophyte %in% "smithora")
zostera<- subset(data,data$macrophyte %in% "zostera")

smithzost <- merge(smith, zostera, by = c("blade_sample_id", "year","date", "site"))
smithzost <- smithzost[,c("site", "date", "blade_sample_id", "macrophyte.x","final_dry_g.x", "macrophyte.y", "final_dry_g.y")]
smithzost$ratio <- smithzost$final_dry_g.x / smithzost$final_dry_g.y
smithzost$date <- as.Date(smithzost$date, format="%d-%b")
#smithzost <- na.omit(smithzost)
#smithzost2 <- subset(smithzost, smithzost$ratio > 0)  

plot(ratio ~ site, data = smithzost)

library(ggplot2)
ggplot( data=smithzost, aes(y=ratio, x=date)) + geom_point() + facet_wrap(~site)
ggplot( data=smithzost, aes(y=ratio, x=date)) + geom_point() + facet_grid(site~.)



