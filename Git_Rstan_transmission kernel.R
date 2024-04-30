# To use STAN for the estimation of transmission kernel parameter -----
# Using data from FMD from Lamphayaklang subdistrict and Boh Ploi district

# Load required package --------
packages <- c("dplyr", "svMisc","igraph", "readxl", "rstan")
lapply(packages, library, character.only = TRUE)

# Load data ---------
LP_distance <- read_excel("LP_distance.xlsx") # between distance between farm in Lamphaya Klang
LP_distance <- LP_distance/1000 # change unit of distance between farms from m to km
BP_distance <- read_excel("BP_distance.xlsx") # between distance between farm in Boh ploi
BP_distance <- BP_distance/1000
kernel_LP <- readRDS("kernel_LP.xlsx") # outbreak data in Lamphaya Klang
kernel_BP <- readRDS("kernel_BP.xlsx") # outbreak data in in Boh ploi

##~~~~~~~~~~~~~~~~~---------------------------------------------
# Manage data to join outbreak data from 2 areas together ------
## For LP-----
sim_data <- kernel_LP
distancematrix <- as.matrix(LP_distance[,-1]) # remove ID row

D = ncol(sim_data) - 2 # only count date
N  = nrow(sim_data) # number of farms

# take the dataframe only with the farm status
sim_data <- sim_data[,-c(1:2)]

# create susceptible matrix
status_susceptible <- as.matrix(sim_data == "S")*1

# create infectious matrix
status_infectious <- as.matrix(sim_data == "I")*1

# create infection event; infection event happen on the last day that it is susceptible
a <- apply(sim_data, 1, function(x) max(which(x %in% c("S"))))
a[a == -Inf] <- ncol(sim_data) # change -Inf (index cases) to the last day

# create infected event matrix
infected_event <- matrix(0, ncol = ncol(status_infectious), nrow = nrow(status_infectious) )
for (i in 1:length(a)){
  infected_event[i, a[i]] <- 1
}
# change the last column to zero; infection event can't happen on the last day
infected_event[, ncol(infected_event)] <- 0

# Create vector of infection event 
new_infect <- c(infected_event)


# create list of force of infection
transmission_event <- list()

for (i in 1:ncol(status_susceptible)){
  
  transmission_event[[i]] <-  (as.matrix(status_susceptible[,i])%*% as.matrix(t(status_infectious[,i])))
  
}

dim(transmission_event[[1]])

# multiply distance to event 
transmission_event <- lapply(transmission_event, function(M) {M * distancematrix})
transmission_event[[1]]* distancematrix
# make data into matrix
transmission_matrix_LP <- do.call(rbind,transmission_event) # every 502 rows represent transmission event of one day

# cbind the event and force of infection
transmission_matrix_LP <- cbind(new_infect, transmission_matrix_LP)
# to check
# which(new_infect == 1 )
# transmission_matrix [3069,1]
nrow(transmission_matrix_LP); ncol (transmission_matrix_LP)
# transpose so distance is in same column, the first row is infection event
transmission_matrix_LP <- t(transmission_matrix_LP)
nrow(transmission_matrix_LP); ncol (transmission_matrix_LP)

# remove column with only zero
transmission_matrix_LP  <- transmission_matrix_LP[, which(colSums(transmission_matrix_LP) > 0)]
nrow(transmission_matrix_LP); ncol (transmission_matrix_LP)

# event
event_LP <- as.numeric(transmission_matrix_LP[1,])

# remove first event row
transmission_matrix_LP <-transmission_matrix_LP[-1,]

# count number of distance != 0 in each column
cut_LP <- colSums(transmission_matrix_LP>0)

# add distance to one vector
distance_kernel_LP <- c(transmission_matrix_LP)

# remove zero
distance_kernel_LP <-distance_kernel_LP [distance_kernel_LP > 0]

# create start and stop for segment
stop_LP <- cumsum(cut_LP)

start_LP <- (cumsum(cut_LP)) - cut_LP+1

## For BP-----
sim_data <- kernel_BP
distancematrix <- as.matrix(BP_distance[,-1]) # remove ID row

D = ncol(sim_data) - 2 # only count date
N  = nrow(sim_data) # number of farms

# take the dataframe only with the farm status
sim_data <- sim_data[,-c(1:2)]

# create susceptible matrix
status_susceptible <- as.matrix(sim_data == "S")*1

# create infectious matrix
status_infectious <- as.matrix(sim_data == "I")*1

# create infection event; infection event happen on the last day that it is susceptible
a <- apply(sim_data, 1, function(x) max(which(x %in% c("S"))))
a[a == -Inf] <- ncol(sim_data) # change -Inf (index cases) to the last day

# create infected event matrix
infected_event <- matrix(0, ncol = ncol(status_infectious), nrow = nrow(status_infectious) )
for (i in 1:length(a)){
  infected_event[i, a[i]] <- 1
}
# change the last column to zero; infection event can't happen on the last day
infected_event[, ncol(infected_event)] <- 0

# Create vector of infection event 
new_infect <- c(infected_event)


# create list of force of infection
transmission_event <- list()

for (i in 1:ncol(status_susceptible)){
  
  transmission_event[[i]] <-  (as.matrix(status_susceptible[,i])%*% as.matrix(t(status_infectious[,i])))
  
}

dim(transmission_event[[1]])

# multiply distance to event 
transmission_event <- lapply(transmission_event, function(M) {M * distancematrix})
transmission_event[[1]]* distancematrix
# make data into matrix
transmission_matrix_BP <- do.call(rbind,transmission_event) # every 502 rows represent transmission event of one day

# cbind the event and force of infection
transmission_matrix_BP <- cbind(new_infect, transmission_matrix_BP)

# transpose so distance is in same column, the first row is infection event
transmission_matrix_BP <- t(transmission_matrix_BP)

# remove column with only zero
transmission_matrix_BP  <- transmission_matrix_BP[, which(colSums(transmission_matrix_BP) > 0)]

# event
event_BP <- as.numeric(transmission_matrix_BP[1,])

# remove first event row
transmission_matrix_BP <-transmission_matrix_BP[-1,]

# count number of distance != 0 in each column
cut_BP <- colSums(transmission_matrix_BP>0)

# add distance to one vector
distance_kernel_BP <- c(transmission_matrix_BP)

# remove zero
distance_kernel_BP <-distance_kernel_BP [distance_kernel_BP > 0]

# create start and stop for segment
stop_BP <- cumsum(cut_BP)
start_BP <- (cumsum(cut_BP)) - cut_BP+1


## Combine LP and BP data  ------
distance_kernel <- c(distance_kernel_LP,distance_kernel_BP)
event <- c(event_LP,event_BP)
cut <- c(cut_LP, cut_BP)
tail(start_LP); tail(stop_LP);tail(cut_LP)
head(start_BP); head(stop_BP);head(cut_BP)
# combine sequence start and stop from LP and BP
a <-start_BP + stop_LP[length(stop_LP)]
b <- stop_BP + stop_LP[length(stop_LP)]
start <- c(start_LP,a)
stop <- c(stop_LP,b)

cut <- c(cut_LP, cut_BP)
K <- sum(cut) # sum cut should be equal to length(distance_kernel)
ID_day <- length(cut)

# Number of farm in both areas
N <- nrow(kernel_LP)+nrow(kernel_BP)


## Run STAN model-----
data <- list(distance_kernel = distance_kernel,event = event, N=N, ID_day =ID_day, start= start, stop = stop, K=K)
require(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

fit <- stan("Git_pilot_transmission_parameter.stan",data = data)
plot (fit)
traceplot(fit)

# check model fit
check_treedepth(fit)
check_energy(fit)
check_divergences(fit)

# extract parameters
k0 <- extract(fit, 'k0')
k0<- unlist(k0, use.names=FALSE)
hist(k0, freq=FALSE)
mean(k0); quantile(k0, 0.025); quantile(k0, 0.975)

r0 <- extract(fit, 'r0')
r0<- unlist(r0, use.names=FALSE)
hist(r0, freq=FALSE)
mean(r0); quantile(r0, 0.025); quantile(r0, 0.975) 

alpha <- extract(fit, 'alpha')
alpha<- unlist(alpha, use.names=FALSE)
hist(alpha, freq=FALSE)
mean(alpha); quantile(alpha, 0.025); quantile(alpha, 0.975) 


