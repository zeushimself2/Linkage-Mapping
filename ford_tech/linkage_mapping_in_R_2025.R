rm(list = ls())

setwd("C:/Users/ZEUS/Desktop/ford_tech")
list.files()

source('custom_functions.R')

# Install and load packages if not already installed
if (!require('googlesheets4')) install.packages('googlesheets4')
if (!require('onemap')) devtools::install_github("augusto-garcia/onemap")
# if (!require('mappoly')) install.packages('mappoly')

# Import BC1 data from google drive
# No authorization to access a shared Google sheet
googlesheets4::gs4_deauth()
# Authorize access to Google Drive
gs4_auth()


# Let's play with some example dataset url to data
bc1_url <- "https://docs.google.com/spreadsheets/d/1gybGwUpItL9f--xIalJFmlrBwT4D3W-CTG5_aX3X4m8/edit?usp=sharing"


# Import clover_data_CRD as a data frame
bc1 <- as.data.frame(googlesheets4::read_sheet(bc1_url,
                                               sheet = 1,
                                               na = "-"))
View(bc1)

# Replace all NA's with a hyphen -- MapMaker requirement
bc1[is.na(bc1)] <- '-'

# Convert backcross data to MapMaker format and save as .txt file
proc_mapmaker(x = bc1,
              type = "f2 backcross",
              geno_start = 2,
              trait_id =  c('PH_cm', 'SD_pct'),
              out_name = "bc_data",
              ext = '.txt')
# The function proc_mapmaker() will create a file called "bc_data.txt" in the working directory

# Import backcross data prepared for MapMaker
dat <- onemap::read_mapmaker(dir = getwd(),
                             file = "bc_data.txt")



dat

# Plot data as a heatmap
plot(dat)


# Set linkage mapping parameters
map_fun <- 'kosambi'

set_map_fun(type = map_fun)


# Find redundant markers
bins <- find_bins(dat, exact = FALSE)
bins

bins_example <- create_data_bins(dat, bins)

# Segregation distortion test for markers
seg_test <- test_segregation(bins_example)
print(seg_test)
plot(seg_test)

# Select distorted markers
select_segreg(seg_test, distorted = TRUE)

# Show the markers numbers without segregation distortion
no_dist <- select_segreg(seg_test, distorted = FALSE, numbers = TRUE)
no_dist

# Show the markers numbers with segregation distortion
dist <- select_segreg(seg_test, distorted = TRUE, numbers = TRUE)
dist

# Suggest LOD for grouping markers
(LOD_sug <- suggest_lod(bins_example))


# Estimate pairwise recombination fractions using default values
twopts_rf <- rf_2pts(input.obj = bins_example,
                     LOD = LOD_sug)

print(twopts_rf, c(1, 2)) # View rf between marker 1 and 2

# Assigning markers to linkage groups
mark_all <- make_seq(twopts_rf, "all")

# Form linkage groups
lg_grps <- group(mark_all)

lg_grps$n.groups # get number of linkage groups

# # Group markers based on clustering
# lg_upgma <- group_upgma(mark_all, expected.groups = 5, inter = FALSE)
# plot(lg_upgma)

lg_maps <- vector(mode = 'list', length = lg_grps$n.groups)

for (i in seq_len(lg_grps$n.groups)) {

  ##create linkage group i
  LG1 <- onemap::make_seq(lg_grps, i)

  # Order linkage groups using different methods
  map.cur <- onemap::ug(LG1, hmm = FALSE)

  lg_maps[[i]]  <- onemap::map(map.cur)

  ##assign the map of the i-th group to the maps.list
  # lg_maps[[i]] <- make_seq(map.cur, "force")


}


draw_map(lg_maps, names = TRUE, cex.mrk = 0.8, cex.grp = .75)

##write maps.list to ".map" file
write_map(lg_maps, "lg_map.txt")

# Draw a pretty map
draw_map2(lg_maps,
          cex.label = 0.5,
          tag = "all",
          main = "Linkage map for BC data",
          group.names = c("LG1", "LG2", "LG3"),
          output = file.path(getwd(), "lg.png"))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ QTL Mapping ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Import linkage map into R as a data frame using the read.table function
linkage_map <- read.table('lg_map.txt', header = FALSE, sep = "")

# # Re-format geno data after linkage mapping to re-order the marker loci
nn <- reorder_markers(geno_dat = bc1,
                      linkage_map = linkage_map,
                      type = "f2 backcross",
                      trait_id =  c('PH_cm', 'SD_pct'),
                      out_name = "bc_data2",
                      ext = '.txt' )


