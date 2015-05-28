#### File to create total Atlanta dataset


# Set output
save.fp <- '~/Dropbox/Dist_exp_sheds/results'
howard <- "/Volumes/data/bios/Howard/For Jenna/"
jenna <- "/Volumes/data/bios/Jenna Krall/SHEDS"

# Load libraries
library(dplyr)
library(reshape2)
library(ncdf4)

# Load data
load(file.path(howard, 'Individual_SHEDS_wID_cleaned.RData'))



# Clean data
dates <- as.numeric(as.Date(dat$date ,format = "%d%b%y"))
dat <- mutate(dat, date = dates)
dat <- select(dat, personID, date, zipcode, pmAmbient)


# Make nc file
dim1 <- ncdim_def("Obs", "obs", 1 : nrow(dat))
ID <- ncvar_def("ID", "", dim1, prec = "integer")
date <- ncvar_def("date", "since1970jan01", dim1, prec = "integer")
zip <- ncvar_def("zip", "", dim1, prec = "integer")
pm <- ncvar_def("pm", "AMBIENTmugm3", dim1, prec = "double")

nc <- nc_create(file.path(jenna, "SHEDS-ncdf.nc"), list(ID, date, zip, pm))
ncvar_put(nc, ID, dat$personID)
ncvar_put(nc, date, dat$date)
ncvar_put(nc, zip, dat$zipcode)
ncvar_put(nc, pm, dat$pmAmbient)
nc_close(nc)
