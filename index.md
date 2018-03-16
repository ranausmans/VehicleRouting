## Localised Vehical Routing with Vendor Based Selection in R

Hello everyone,


For past few days, I developed a strong interesting in building a basic implementation of Vehical routing in R. The problem can be approach through many ways as there is no single handed solution to it. Vehical Routing is extension of famous algorithmic problem called TSP or Traveling Salesman Problem. It is to optimize the route for a traveling salesman to reach maximum locations in less time.

I have worked on extension of TSP for Vehical Routing with addition to selection of Vendor/Truck Depot for deliveries closest to the Depot. 

### Scenario

Consider there are mulitple depots and multiple customers. In the image shown below, I-10 and F-11 are addresses of Depots and rest of addresses belong to Customers.

![img1](https://i.imgur.com/hq5kWZo.png)

Currently, in a nutshell, the addresses are scattered.

### Finding Geocodes for addresses.

The code is self explanatory. It takes in the csv as shown in above image and finds the geo location using Google API.

```

#load ggmap
library(ggmap)

origAddress <- read.csv2("G:/funproj/docs/adds.csv", header=FALSE, sep="")

# Initialize the data frame
geocoded <- data.frame(stringsAsFactors = FALSE)

# Loop through the addresses to get the latitude and longitude of each address and add it to the
# origAddress data frame in new columns lat and lon
for(i in 1:nrow(origAddress))
{
  # Print("Working...")
  result <- geocode(origAddress$addresses[i], output = "latlona", source = "google")
  origAddress$lon[i] <- as.numeric(result[1])
  origAddress$lat[i] <- as.numeric(result[2])
  origAddress$geoAddress[i] <- as.character(result[3])
}
# Write a CSV file containing origAddress to the working directory
write.csv(origAddress, "geocoded.csv", row.names=FALSE)
```

This will annotate the latitudes and longitudes, the output would look like this.

![img2](https://i.imgur.com/cr25Npp.png)

The next thing is to cluster the data by number of Vendors in such a way that the closest customers to Depot are aligned.

At this point, the data is scattered.

![img3](https://i.imgur.com/s6r97mY.png)

But after you run the following code.

```

library(TSP)
library(stringi)

location <- read.csv("geocodes.csv", header=FALSE, sep="")

# Compute locations in polar coordinates (to work better for euclidean distances)
location$LONG_RAD <- location$LONGITUDE * (2 * pi)/360
location$LAT_RAD <- (location$LATITUDE * 2) * (2 * pi)/360
R <- (6378 + 6356)/2
location$X = R * cos(location$LAT_RAD) * cos(location$LONG_RAD)
location$Y = R * cos(location$LAT_RAD) * sin(location$LONG_RAD)
location$Z = R * sin(location$LAT_RAD)

# Group cusomer addresses into n clusters, where n is the number of representative addresses
n_cluster <- length(location$TYPE[location$TYPE == "VENDOR"])
group1 <- kmeans(x=location[location$TYPE=="CUSTOMER", colnames(location) %in% c("X", "Y", "Z")], centers=n_cluster)


closest_cluster <- function(x, c) {
  cluster_dist <- apply(c, 1, function(y) sqrt(sum((x-y)^2)))
  return(which.min(cluster_dist)[1])
}

# Assign each representative to its closest customer cluster
# Not an optimal solution, but it is a good approximation
group2 <- NULL
centers <- group1$centers
addresses <- location[location$TYPE=="VENDOR", colnames(location) %in% c("X", "Y", "Z")]
for(i in 1:(nrow(addresses)-1)) {
  address <- addresses[i, ]
  closest <- closest_cluster(address, centers[, colnames(centers) %in% c("X", "Y", "Z")])
  group2 <- c(group2, as.integer(names(closest)))
  centers <- centers[-closest,]
}
group2 <- c(group2, setdiff(1:nrow(group1$centers), group2))

location$GROUP <- c(group1$cluster, group2)
location$ID <- paste0("C", 1:nrow(location))
location$ID[location$TYPE=="VENDOR"] <- paste0("R", 1:n_cluster)

write.csv(location, locs.csv)
```

At this point, the Vendors are grouped with closest customers, so something of a sort. I run K-mean to assign two groups based on Vendors and closest distance between Vendor and Customer Locations based on Geocodes. In the process, I convert geocodes into Radians to work better with Euclidean distances, a primary base for K-mean.

![img4](https://i.imgur.com/EUDJ4fx.png)

The third script is to run the schedule using the Farthest Insertion Algorithm. 

```
library(TSP)
library(stringi)
library(dplyr)


location <- read.csv("locs.csv", header=FALSE, sep="")

# Compute optimal path (shortest euclidean distance) for each customer address cluster
paths <- list()
distances <- list()
n_cluster <- length(location$TYPE[location$TYPE == "VENDOR"])
for(i in 1:n_cluster) {
  x <- dist(location[location$GROUP==i, colnames(location) %in% c("X","Y","Z")], method="euclidean", diag=T, upper=T)
  x1 <- as.matrix(x)
  rownames(x1) <- location[location$GROUP==i, "ID"]
  colnames(x1) <- location[location$GROUP==i, "ID"]
  tsp <- as.TSP(x1)
  atsp <- as.ATSP(tsp)
  v <- which(stri_startswith_fixed(labels(tsp), "R"))
  atsp[, v] <- 0
  path <- solve_TSP(atsp, method="farthest_insertion")
  path <- cut_tour(path, v, exclude_cut=F)
  paths[[i]] <- path
  distances[[i]] <- x1
}

# Organize the data for each computed path to present the information in terms of start location, end location, sequence number, and distance
data_path <- data.frame()
for(i in 1:length(paths)) {
  for(j in 1:(length(paths[[i]])-1)) {
    d <- distances[[i]]
    from <- labels(paths[[i]][j])
    to <- labels(paths[[i]][j+1])
    type <- "CUSTOMER"
    if(stri_startswith_fixed(from, "R")) {
      type <- "VENDOR"
    }
    lat_from <- location$LATITUDE[location$ID == from]
    lon_from <- location$LONGITUDE[location$ID == from]
    lat_to <- location$LATITUDE[location$ID == to]
    lon_to <- location$LONGITUDE[location$ID == to]
    data_path <- rbind(data_path, data.frame(GROUP=i, FROM=from, TO=to, SEQUENCE=j, DISTANCE=d[from, to], TYPE=type,
                                             LAT_FROM=lat_from, LON_FROM=lon_from, LAT_TO=lat_to, LON_TO=lon_to))
  }
}

```
The output is like this.

![img5](https://i.imgur.com/YJfjz5C.png)

It has scheduled the best optimal path and if someone who lives in Islamabad, Pakistan, he would know that this would be the best way to move forward. 

I don't say this is perfect implementation as Vehical Routing is far more serious and deeper problem, but this is good enough for someone who wants to get started with VRP. The improvements can be made using any evolutionary algorithm such as Genetic Algorithm or Ant Colony but it is to notice that all evolutionary algorithms are NP-Hard. 

I would be happy if someone can create Shiny Application with this basic implementation. 

*Rana Muhammad Usman is a data scientist based in Islamabad, Pakistan who works with Persontyle Group. He has worked on various problems in Data Science such as in Oil Industry and Web User Behaviorial Analysis. He likes to run and aspires to participate in 10k Marathon. He can be reached out at usmanashrafrana@gmail.com*
