
setwd("~/Desktop/stat202a/proj")
rm(ls()) # clear session
options(scipen = 10) # avoid scientific notation
pacman::p_load(sf, tidycensus, dplyr, ggplot2, stringr)
CENSUSKEY <- "e5bda6570015f98cd0bba050c9d21d613d7d4f55"

# From https://data.lacity.org/Public-Safety/Crime-Data-from-2020-to-Present/2nrs-mtv8
la <- read.csv("data/Crime_Data_from_2020_to_Present.csv")
la <- la %>% filter(LAT != 0, LON != 0)

# from https://data.lacity.org/api/views/2nrs-mtv8/files/787064ac-a2f2-4cf3-a474-45036b7da937?download=true&filename=UCR-
# COMPSTAT062618.pdf
carcodes <- c(510,520,330,331,410,420,421) # includes vehicle thefts
car2codes <- c(330,331,410,420,421) # only break ins
bikecodes <- c(480,485) 
homicidecodes <- c(110,113)

bikepts <- la[la$Crm.Cd.1%in% bikecodes,] %>% st_as_sf(., coords = c("LON","LAT"), crs=4326)
carpts <- la[la$Crm.Cd.1%in% carcodes,] %>% st_as_sf(., coords = c("LON","LAT"), crs=4326)
carpts2 <- la[la$Crm.Cd.1%in% car2codes,] %>% st_as_sf(., coords = c("LON","LAT"), crs=4326)

murderpts <- la[la$Crm.Cd.1%in% homicidecodes,] %>% st_as_sf(., coords = c("LON","LAT"), crs=4326)
dim(carpts)/ dim(la) - dim(carpts2)/ dim(la)
dim(bikepts)/ dim(la)


# Plot frequency of crimes
la <- la %>% mutate(crimecat = ifelse(Crm.Cd.1 %in% carcodes, "MTV/BFTV",
                                      ifelse(Crm.Cd.1 %in% c(121,122,815,820,821), "rape",
                                             ifelse(Crm.Cd.1 %in% c(310,320,210,220), "robbery/burglary",
                                                    ifelse(Crm.Cd.1 %in% c(235,236,250,251,761,926,626,627,647,763,928,930), "domestic violence",
                                                           ifelse(Crm.Cd.1 %in% c(230,231,435:437,622:625), "assault", 
                                                                  ifelse(Crm.Cd.1 %in% c(350:353,450:453, 341:345, 440:445, 470:475), "theft", "other")))))),
                    crimecat = factor(crimecat, levels=names(sort(table(crimecat), decreasing=TRUE))) )
ggplot() + geom_bar(data = la %>% filter(crimecat !="other"), aes(crimecat)) + xlab("")

# I only mutate subsets of the full crime data bc og file is too large
addTimeCols <- function(df) {
  return(df %>% mutate(date = str_extract(DATE.OCC, "\\S+"),
                           date = as.Date(date, "%m/%d/%Y"),
                           month = lubridate::month(date),
                           season = ifelse(month %in% 3:5, "spring",
                                           ifelse(month %in% 6:8, "summer",
                                                  ifelse(month %in% 9:11, "fall", "winter"))),
                           hour = ifelse(nchar(TIME.OCC) == 2, 0,
                                         ifelse(nchar(TIME.OCC) < 4, substr(TIME.OCC,1, 1), substr(TIME.OCC, 1, 2))),
                           hour = as.numeric(hour),
                           timeofday = ifelse(hour %in% c(0:5), "midnight",
                                              ifelse(hour %in% 6:11, "morning",
                                                     ifelse(hour %in% 12:17, "afternoon", "night")))))
}
carpts <- addTimeCols(carpts)
murderpts <- addTimeCols(murderpts)

# Temporal trends
# mtv
library(gridExtra)
p1 <- ggplot() + geom_bar(data = carpts, aes(hour)) + scale_x_continuous(breaks=0:23) + ggtitle("Number of car thefts by hour of day")
p2 <- ggplot() + geom_bar(data = carpts, aes(month)) + scale_x_continuous(breaks=1:12)+ ggtitle("Number of car thefts by month")
grid.arrange(p1, p2, ncol=2)

#murder
ggplot() + geom_bar(data = murderpts, aes(hour)) + scale_x_continuous(breaks=0:23)
ggplot() + geom_bar(data = murderpts, aes(month)) + scale_x_continuous(breaks=1:12)


##### Predicting count at neighborhood level

# if you're judging me for how much dplyr I use, I don't care I think it's easier to read
# ACS
vars <- c("B01003_001", "B25046_001","B19326_001","B25008_002", "B25008_003", "B25077_001" )
names(vars) <- c("TotalPop", "NumVehicles", "MedInc", "Owners", "Renters", "MedHomeVal")
acs <- get_acs(geography = "tract", variables = vars,
               state = "CA", county = "Los Angeles", 
               year = 2019, geometry = T, cache_table = TRUE,
               key = CENSUSKEY,  survey = "acs5") %>% select(-NAME)
acs <- st_transform(acs, crs = 4326)

# Interpolate missing data for tracts by using average among neighboring tracts
centers <- st_centroid(acs[!duplicated(acs$GEOID),], byid = T)
centers <- centers %>% cbind(st_coordinates(centers))

# Some tracts have missing data, use average of 6 nearest neighbors to interpolate its value
# I use distance between tract centroids to define nearest neighbors
for (var in unique(acs$variable)) {
  rowsubset <- acs$variable==var 
  need2fix <- acs[rowsubset & is.na(acs$estimate),]$GEOID
  for(geoid in need2fix) {
    tractX <- centers[centers$GEOID == geoid,]$X
    tractY <- centers[centers$GEOID == geoid,]$Y
    dists <- (centers$X - tractX)^2 + (centers$Y - tractY)^2
    # 6 nearest neighbors + 1 including the tract itself, its NA value is dropped from the mean
    nearest <- centers[head(order(dists), 6+1),]$GEOID
    estimate <- mean(acs[(acs$GEOID %in% nearest)&rowsubset, ]$estimate, na.rm = T)
    acs[(acs$GEOID==geoid) & rowsubset, "estimate"] <- estimate
  }
}

dat <- reshape(data.frame(acs) %>% select(-geometry), idvar = "GEOID", timevar = "variable", direction = "wide")

getCounts <- function(codes, lab) {
  sub <- la[la$Crm.Cd.1%in% codes,]
  pts <- sub[,c("Date.Rptd","Crm.Cd.Desc","LAT","LON")] %>% st_as_sf(coords = c("LON","LAT"), crs = 4326)
  countbytract <- (st_join(pts, acs) %>% as.data.frame() %>% group_by(GEOID) %>% summarise(lab0 = n()))
  names(countbytract)[names(countbytract) == "lab0"] <- lab
  return(countbytract)
}

dat <-( dat %>% merge(., getCounts(carcodes, "numcarthefts"), by = "GEOID", all.y = T)
        %>% merge(., getCounts(bikecodes, "numbikethefts"), by = "GEOID", all.x = T)
        %>% merge(., getCounts(homicidecodes, "numhomicides"), by = "GEOID", all.x = T))

  
# Valid zipcodes
zipcodes <- st_read("data/lA_zipcodes/", crs = 4326)
usedzipcodes <-(zipcodes %>% st_join(carpts, .) 
                %>% group_by(zipcode) %>% summarise(numcarthefts = n())
                %>% filter(numcarthefts > 0) %>% pull(zipcode))
usedzipcodes <- as.character(usedzipcodes)
# areas with 0 thefts reported weren't enforced by LAPD

## Businesses per tract
# From https://data.lacity.org/Administration-Finance/Listing-of-Active-Businesses/6rrh-rzua/data
bus <- read.csv("data/Listing_of_Active_Businesses.csv")
dim(bus)
bus <- (bus %>% select(ZIP.CODE, BUSINESS.NAME, NAICS, PRIMARY.NAICS.DESCRIPTION, COUNCIL.DISTRICT, LOCATION)
        %>% mutate(ZIP.CODE = gsub("-\\d+", "", ZIP.CODE))
        %>% mutate(LAT = as.numeric(str_extract(LOCATION, "(\\d{2}.\\d+)")),
                   LON = as.numeric(str_extract(LOCATION, "(-\\d{2,3}.\\d+)")))
        %>% filter(!is.na(LON), !is.na(LAT)))
dim(bus)
buspts <- bus[,c("LAT","LON")] %>% st_as_sf(coords = c("LON","LAT"), crs = 4326)
countbytract <- st_join(buspts, acs) %>% as.data.frame() %>% group_by(GEOID) %>% summarise("BusinessCount" = n())
dat <- merge(dat, countbytract, by = "GEOID", all.x = T)

# **** Final dat ****
# currently numbikethefts/numhomicides=NA when there were none reported in that tract
dat <- dat %>% mutate(across(where(is.numeric), function(v)  ifelse(is.na(v), 0, v)))
names(dat) <- gsub("estimate.", "", names(dat))
dat$RenterProp <- dat$Renters / dat$TotalPop

#### GAM
library(mgcv)
g <- gam(numcarthefts ~ s(BusinessCount) + s(TotalPop) + s(MedInc) #+ s(estimate.Owners) 
         + s(RenterProp) + s(MedHomeVal), data = dat)
summary(g)
plot(g, residuals = T, cex = 3, shift = coef(g)[1], ylab = "")

gmurder <- gam(numhomicides ~ s(BusinessCount) + s(TotalPop) + s(MedInc) #+ s(estimate.Owners) 
         + s(RenterProp) + s(MedHomeVal), data = dat)
summary(gmurder)
plot(gmurder, residuals = T, cex = 3, shift = coef(g)[1], ylab = "")


##### Kernel smoothing 
dyn.load("~/Desktop/stat202a/proj/proj.so")
kde2dR <- function(x, grd) {
  bw <- apply(x, 2, function(v) bw.nrd(v[!is.na(v)])) 
  bw <- sqrt(sum(bw^2))
  n <- nrow(x) %>% as.integer()
  p <- ncol(x) %>% as.integer()
  m <- nrow(grd) %>% as.integer()
  out <- .C("mykde2d", n, p, m, bw, t(as.matrix(grd)),t(as.matrix(x)), dens = double(m))
  return(out$dens)
}

# I want reference labels bc I don't know LA by shape so i'll label a sample of neighborhoods... aka me doing the most
keep <- c("Westwood", "Boyle Heights", "Venice", "Pacific Palisades",
          "Hollywood", "Pacoima", "Koreatown", "Florence", "Tujunga",
          "Valley Village", "Echo Park", "Reseda", "Granada Hills",
          "Highland Park", "Mar Vista", "Wilmington", "Downtown",
          "San Pedro", "Sylmar", "Panorama City", "Palms", "Vermont Square",
          "Uniersity Park", "Century City", "Mid-City", "Sherman Oaks")
neighLabs <- (st_read("data/laneighs/") %>% st_centroid())
neighLabsF <- (data.frame(st_coordinates(neighLabs), name = neighLabs$name)
               %>% filter(name %in% keep))

gridpts <- zipcodes %>% filter(zipcode %in% usedzipcodes) %>% st_make_grid(n = 100, what = "centers")
kernels <- data.frame(st_coordinates(carpts)) 
kernelsbike <- data.frame(st_coordinates(bikepts))
gridpts <- data.frame(st_coordinates(gridpts))
# density based on proximity to other thefts
gridpts$density <- kde2dR(kernels, gridpts %>% select(X, Y)) # density at a point location
gridpts$bikedensity <- kde2dR(kernelsbike, gridpts %>% select(X, Y))


# Map
(ggplot() 
  + geom_point(data = gridpts, aes(x=X, y=Y, color = bikedensity))
  + scale_color_continuous(low = "gray", high = "red")
  + geom_text(data = neighLabsF, aes(x=X, y=Y, label=name), size = 3)
  + xlab("") + ylab("")
)

gridpts %>% arrange(-bikedensity) %>% head(5)

# Differences between seasons?
gridpts2 <- zipcodes %>% st_make_grid(n = 100, what = "centers")
gridpts2 <- data.frame(st_coordinates(gridpts2))

for (s in unique(carpts$season)) {
  kerns <- data.frame(st_coordinates(carpts %>% filter(season == s))) 
  gridpts2[,s] <- kde2dR(kerns, gridpts2 %>% select(X, Y))
}

gridpts2 <- gridpts2 %>% reshape2::melt(id.vars = c("X","Y"), variable.name = "season", value.name= "density")

(ggplot() 
  + geom_point(data = gridpts2, aes(x=X, y=Y, color = density))
  + scale_color_continuous(low = "gray", high = "red")
  + xlab("") + ylab("") 
  + facet_wrap(~ season)
)

# contrast to murder...
gridpts2 <- zipcodes %>% st_make_grid(n = 100, what = "centers")
gridpts2 <- data.frame(st_coordinates(gridpts2))
for (s in unique(kernels2$season)) {
  kerns <- data.frame(st_coordinates(murderpts %>% filter(season == s))) 
  gridpts2[,s] <- kde2dR(kerns, gridpts2 %>% select(X, Y))
}

gridpts2 <- gridpts2 %>% reshape2::melt(id.vars = c("X","Y"), variable.name = "season", value.name= "density")

(ggplot() 
  + geom_point(data = gridpts2, aes(x=X, y=Y, color = density))
  + scale_color_continuous(low = "gray", high = "red")
  + xlab("") + ylab("") 
  + facet_wrap(~ season)
)

### business data eda
bigcats <- bus %>% group_by(PRIMARY.NAICS.DESCRIPTION) %>% summarise(n = n()) %>% arrange(-n) 
bigcats$PRIMARY.NAICS.DESCRIPTION %>% head(25)
catty <- bigcats$PRIMARY.NAICS.DESCRIPTION[2]
sub <- bus %>% filter(PRIMARY.NAICS.DESCRIPTION == catty,!is.na(LON),!is.na(LAT))
head(sub$BUSINESS.NAME)
gridpts$busdensity <- kde2dR(sub %>% select(LON,LAT), testpts%>% select(X,Y))





