if(!require("osmdata")) install.packages("osmdata")
if(!require("tidyverse")) install.packages("tidyverse")
if(!require("sf")) install.packages("sf")
if(!require("ggmap")) install.packages("ggmap")

library(tidyverse)
library(osmdata)
library(sf)
library(ggmap)

# https://dominicroye.github.io/en/2018/accessing-openstreetmap-data-with-r/


# https://wiki.openstreetmap.org/wiki/Map_Features
available_tags("highway")
available_tags("waterway")

# https://www.rdocumentation.org/packages/osmdata/versions/0.1.1/topics/getbb
getbb("Wien")
city <- getbb("Wien")

# https://rdrr.io/cran/osmdata/man/opq.html
# https://github.com/ropensci/osmdata
# https://www.rdocumentation.org/packages/osmdata/versions/0.1.1/topics/osmdata_sf
streets <- opq(city) %>%
  add_osm_feature(key = "highway",
                  value = c("motorway", "primary", "secondary", "tertiary")) %>%
  osmdata_sf()
streets

small_streets <- opq(city) %>%
  add_osm_feature(key = "highway",
                  value = c("residential", "living_street", "unclassified", "service", "footway")) %>%
  osmdata_sf()

river <- opq(city) %>%
  add_osm_feature(key = "waterway",
                  value = "river") %>%
  osmdata_sf()

# https://ggplot2.tidyverse.org/reference/ggsf.html
p1 <- ggplot() +
  geom_sf(data = streets$os_lines, inherit.aes = FALSE, 
          color = "red", size = .4, alpha = .8) +
  geom_sf(data = small_streets$osm_lines, inherit.aes = FALSE,
          color = "black", size = .2, alpha = .6) +
  geom_sf(data = river$osm_lines, inherit.aes = FALSE,
          color = "blue", size = .2, alpha = .5) +
  coord_sf(xlim = c(city[1, 1], city[1, 2]), 
           ylim = c(city[2, 1], city[2, 2]),
           expand = FALSE)
p1

p2 <- ggplot() +
  geom_sf(data = streets$os_lines, inherit.aes = FALSE, 
          color = "red", size = .4, alpha = .8) +
  geom_sf(data = small_streets$osm_lines, inherit.aes = FALSE,
          color = "gold", size = .2, alpha = .6) +
  geom_sf(data = river$osm_lines, inherit.aes = FALSE,
          color = "blue", size = .5, alpha = .5) +
  coord_sf(xlim = c(city[1, 1], city[1, 2]), 
           ylim = c(city[2, 1], city[2, 2]),
           expand = FALSE) +
  theme_void() +
  theme(plot.background = element_rect(fill="grey10"))
p2

backmap <- get_map(city, maptype = "toner-background")
p3 <- ggmap(backmap)
p3

cine <- opq(city) %>%
  add_osm_feature("amenity", "cinema") %>%
  osmdata_sf()

p3 +
  geom_sf(data = cine$osm_points, inherit.aes = FALSE, 
          color = "red", size = 2, alpha = .4) +
  labs(x="", y="")
