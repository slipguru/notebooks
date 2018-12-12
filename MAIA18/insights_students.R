library("tidyverse") # load data manipulation library
library("downloader") # load data download utility library

# Move to the desired working directory
working_path <- '..'
setwd(working_path)

# Download data
url <- "https://raw.githubusercontent.com/genomicsclass/dagdata/master/inst/extdata/msleep_ggplot2.csv"
if (!file.exists(basename(url))) download(url, basename(url))

# Load data into a tibble (named msleep) and inspect its content


# Inspect variables description
variables <- tribble(
  ~column_name,	~description, ~um,
  #--|--|----
  "name", "common name", NA,
  "genus",	"taxonomic rank", NA,
  "vore",	"carnivore, omnivore or herbivore?", NA,
  "order",	"taxonomic rank", NA,
  "conservation",	"the conservation status of the mammal", NA,
  "sleep_total",	"total amount of sleep", "hours",
  "sleep_rem",	"rem sleep", "hours",
  "sleep_cycle",	"length of sleep cycle", "hours",
  "awake",	"amount of time spent awake", "hours",
  "brainwt",	"brain weight", "kilograms",
  "bodywt",	"body weight", "kilograms"
  )
print(variables)

# From now on, try to use pipes if possible...
# Inspect all the unique mammals genera


# Keep only the mammals with brain heavier than 0.1 kg and show them sorted by body weight


# Save name, brain and body weight of the mammals above in a new tibble called msmart


# Make scatterplot of brain weight (horizontal) and body weight (vertical) of such mammals
# (take care of naming axes)


# Who are the outliers? Drop them and re-create the plot


# Do carnivores sleep more than herbivores?
# [1] Drop mammals with vore = NA and with less than 10 samples per vore


# [2] Inspect sleep_total distribution


# [3] Run hypothesis test (critical value = 5%)


# replace `???` with the correct R operator
if (t_test$p.value ??? alpha){
  print(paste('Carnivores sleep more than herbivores (alpha =', as.character(alpha),')'))
} else {
  print(paste('Carnivores do NOT sleep more than herbivores (alpha =', as.character(alpha),')'))
}


  
