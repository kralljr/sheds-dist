load('/Volumes/data/bios/Howard/For Jenna/Individual_SHEDS.RData')

save.fp <- '~/Dropbox/Dist_exp_sheds/results'
library(dplyr)
library(reshape2)


dates <- as.Date(dat[, 1] ,format = "%d%b%y")



dat2 <- select(dat, date, zipcode, pmAmbient)
dat2 <- mutate(dat2, date = dates)


#unique dates
un_date <- unique(dates)



# Determine whether normal
pvals <- list()
errs <- c(0, 0)

for(i in 1 : length(un_date)) {
	perc <- round(i / length(un_date), 2)
	print(c(i, perc))
	pvals[[i]] <- vector()
	pmdat1 <- filter(dat2, date == un_date[i])
	un_zip <- unique(pmdat1$zipcode)
	for(j in 1 : length(un_zip)) {
		
		pmdat <- filter(pmdat1, zipcode == un_zip[j])
		temp <- try(shapiro.test(pmdat[, 3])$p.val)
		if(class(temp) == "try-error") {
			errs1 <- c(un_date[i], un_zip[j])
			errs <- rbind(errs, errs1)
			pvals[[i]][j] <- NA
		}else {
			pvals[[i]][j] <- temp
		}
		
	}
	
}


names(pvals) <- un_date
mpvals <- melt(pvals)
colnames(mpvals) <- c("SWpval", "Date")

write.csv(mpvals, row.names = F, file = file.path(save.fp, "SWpvals.csv"))

#reject in all cases??
length(which(mpvals[, 1] < 0.05)) / nrow(mpvals)



pdf(file.path(save.fp, "hist_pm.pdf")) 
par(mfrow = c(3, 3))
for(i in 1 : 9) {
	sam_date <- sample(un_date, 1)
	dat3 <- dat2 %>% filter(date == sam_date)
	sam_zip <- sample(dat3$zipcode, 1)

	pm <- dat3 %>% filter(zipcode == sam_zip) %>% select(pmAmbient)
	name1 <- paste0("Date = ",sam_date, ", Zip = ", sam_zip)
	hist(pm[, 1], main = name1, xlab = "Ambient PM")
}
graphics.off()
