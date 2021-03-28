# Set working directory
setwd("~/Google Drive/HKU Drive/Google Drive/School/2020-21/FYP/D614G/Example_D614G_for_Kay/Data Manipulation")

# Library
library(dplyr)

# Read file
total = read.csv("owid-covid-data.csv")

countryListConfirmed = c("Australia", "Belgium", "Denmark", "Iceland", "India", 
                         "Netherlands", "Spain", "Portugal", "United Kingdom", "United States")

for (country in countryListConfirmed) {
  tmp = total %>% 
    filter(location == country) %>% 
    subset(select = c(date, new_deaths))
  tmp$date = as.Date(tmp$date)
  tmp$day_num = as.numeric(tmp$date - as.Date("2019-12-22") + 1)
  tmp$week_num = tmp$day_num %/% 7 + 1
  if (country == "Spain") {
    n = which(tmp$date == "2020-05-25")
    tmp$new_deaths[n] = (tmp$new_deaths[n-1] + tmp$new_deaths[n+1])/2
  }
  write.csv(tmp, file = paste("../0. data/confirmed_by_country/owid_death_", country, ".csv", sep = ""),
            row.names = FALSE)
}
