####
# File to coarsen SHEDS data

library(ncdf4)
library(dplyr)
library(reshape2)
library(bigtabulate)

# Set output
save.fp <- '~/Dropbox/Dist_exp_sheds/results'
howard <- "/Volumes/data/bios/Howard/For Jenna/"
jenna <- "/Volumes/data/bios/Jenna Krall/SHEDS"

data <- nc_open(file.path(jenna, "SHEDS-ncdf.nc"))
pm <- ncvar_get(data, "pm")
ID <- ncvar_get(data, "ID")
date <- ncvar_get(data, "date")
zip <- ncvar_get(data, "zip")

nc_close(data)

dat <- data.frame(ID, date, zip, pm)
dat <- apply(dat, 2, as.numeric)
dat <- data.frame(dat)
dates <- as.Date(dat$date, origin = "1970-01-01")
dat <- mutate(dat, date = dates)

pop <- read.csv("~/Dropbox/SHEDS/data/Pop_ZCTA_2000.csv", stringsAsFactors = F)
zip <- strsplit(pop[, 1], " ") %>% sapply(., function(x) x[2]) %>% as.numeric()


# Howard's way:
# Largest zip code
zip1 <- zip[which.max(pop[, 2])]
# How many to sample in each zip
pop2 <- ceiling(pop[, 2] / pop[which.max(pop[, 2]), 2] * 100)

# Create matrix to sample from
selpop <- data.frame(zip, pop2 )
selpop <- selpop[complete.cases(selpop), ]


unzip <- unique(dat$zip)




# zipfun <- function(dat, pop2) {
	# newdat <- dat[1, ]
	
	
	# # For each zip code
	# for(i in 1 : length(unzip)) {
		# print(i)
		# # Restrict to data for that zip
		# dat2 <- dplyr::filter(dat, zip == unzip[i])
		# # Number to sample for each day
		# nS <- selpop[selpop[, 1] == unzip[i], 2]
		# # Number of days
		# und <- unique(dat2$date)
		
		# # For each day
		# for(j in 1 : length(und)) {
			# dat3 <- filter(dat2, date == und[j])
			# samps <- sample(seq(1, nrow(dat3)), nS, replace = F)
			# newdat <- rbind(newdat, dat3[samps, ])
		# }
	# }
	# newdat[-1, ]
# }

# newdat <- zipfun(dat, pop2)
# write.csv(newdat, file.path(jenna, "SHEDS_sampled.csv"), row.names = F)


