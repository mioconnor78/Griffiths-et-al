x <- c("ggmap", "rgdal", "rgeos", "maptools", "dplyr", "tidyr", "tmap")
lapply(x, library, character.only = TRUE) 

getwd()
setwd("C:/Users/Gwen/Desktop/Rpractice/mapping_in_R")
library(rgdal)

lnd <- readOGR(dsn = "data", layer = "london_sport")

seagrass <- readOGR(dsn = "data", layer = "points2")

head(lnd@data, n = 2)

head(lnd@data, n = 2)

head(lnd@polygons[[1]]@Polygons[[1]]@coords, 3)

plot(lnd@polygons[[1]]@Polygons[[1]]@coords)

plot(lnd)

plot(seagrass)

sel <- lnd$Partic_Per > 20 & lnd$Partic_Per < 25
plot(lnd[sel, ]) # output not shown here
head(sel) # test output of previous selection (not shown)

plot(lnd, col = "lightgrey") # plot the london_sport object
sel <- lnd$Partic_Per > 25
plot(lnd[ sel, ], col = "turquoise", add = TRUE)

library(leaflet)

m <- leaflet() %>%
  addTiles() %>%  # Add default OpenStreetMap map tiles
  addMarkers(lng=-128.1170464, lat=51.67250839	, popup="The birthplace of R")
m

##Gwen's actual work
sites <- read.csv("choked_sites.csv")
getwd()
df = data.frame(sites)
leaflet(df)  %>% addCircles()

choked <- leaflet(df) %>%
  addTiles() %>%
  addCircleMarkers(radius = ~scale, fill = TRUE, color = "red") %>%
  addLabelOnlyMarkers(~Long, ~Lat, label =  ~as.character(Name), labelOptions = labelOptions(noHide = T, direction = 'top', textOnly = T))

head(sites)
choked
print(choked)

head(sites)

library(DiagrammeR)
?DiagrammeR

DiagrammeR("graph TB;A(Rounded)-->B[Squared];B-->A(Rounded);
           style A fill:#E5E25F;  style B fill:#87AB51;"
)

