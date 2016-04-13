### File to check individual files



# Set directories
setwd("~/Documents/SHEDS")


# Get libraries
library(dplyr)



### Get population counts for each ZCTA
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

######################################






####### Get files of data
lf <- list.files()
# for PM
#lf <- lf[which(substr(lf, 1, 4) == "file")]
cn1 <- c("ID", "date", "zip", "pm")

# for EC
lf <- lf[which(substr(lf, 1, 7) == "EC_file")]
cn1 <- c("ID", "date", "zip", "EC")
######################################



######## Load old data to check length of zips
# load("Individual_SHEDS_wID_cleaned.RData")
# unzip <- unique(dat$zipcode)

#length(unzip) == length(lf)
########





########### Sample up to Atlanta
set.seed(15630)

nrows <- vector()
#zips <- as.numeric(substr(lf, 6, 10))
zips <- as.numeric(substr(lf, 9, 13))

weird1 <- 0
undates1 <- vector()
keeps <- c(0, 0, 0, 0)
pop1 <- vector()


for(i in 1 : length(lf)) {
	k <- round(i / length(lf) * 100, 1)
	
	#print(c(i, k))
	
	# Read in data and create column names
	dat <- read.csv(lf[i], header = F)
	colnames(dat) <- cn1
	
	# # Fix date
	dat$date <- as.numeric(as.Date(dat$date, format = "%d%b%y"))
	
	# Find corresponding population count
	pop <- selpop[selpop[, 1] == zips[i], 2]

	
	# Save nrow (to check against original data)
	nrows[i] <- nrow(dat)
	
	
 	undates <- unique(dat$date)
 	undates1[i] <- length(undates)
 	
	dat <- dat[order(dat$date), ]
	
	# check divide properly
	if((nrow(dat) %% length(undates)) != 0) {
		weird1 <- c(weird1, i)
	}
	
	if(length(unique(dat[, 1])) != (nrow(dat) / length(undates))) {
		print(c(zips[i], i))
	}
	
	s1 <- 1
	
	if(length(pop) != 0) {
		pop1[i] <- pop
	
		#Divide the data by date and randomly sample pop rows
		s1 <- lapply(split(dat, dat[, "date"]), function(x) {
			# randomly shuffle the rows
			x <- x[sample(nrow(x), pop), ]
			rownames(x) <- paste0(seq(1, nrow(x)), x[, 2])
			x
		}) %>% bind_rows %>% as.matrix
				
		keeps <- rbind(keeps, s1)		
	
	}
	rm(dat, s1)
	gc()
	
}

keeps <- keeps[-1, ]
colnames(keeps) <- cn1
keeps <- data.frame(keeps)
keeps[, "date"] <- as.Date(keeps[, "date"], origin = "1970-01-01")
#write.csv(keeps, file = "SHEDS_pm_atl.csv", row.names = F)

write.csv(keeps, file = "SHEDS_ec_atl.csv", row.names = F)


#keeps <- read.csv("SHEDS_pm_atl.csv", stringsAsFactors = F)

# length in all files equals length of total data
