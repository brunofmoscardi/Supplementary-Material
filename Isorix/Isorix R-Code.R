# Load required packages
library(IsoriX)
require(terra)

# Set number of CPUs for IsoriX parallel processing
options_IsoriX(Ncpu = 4)

# Read and prepare water isotope data
rawWater <- read.csv("Water_Nqn_2022F2_sMzaNNqn.csv", h=TRUE)
rawWater$source_ID <- as.factor(rawWater$source_ID)
colnames(rawWater) <- c("source_ID", "lat", "long", "elev", "year", "month", "source_value")
rawWater

# Aggregate and filter water data for study region
rawWateragg <- prepsources(data = rawWater,
                           long_min = -72.012 , long_max = -67.380,
                           lat_min = -39.750, lat_max = -33.919)
rawWateragg

# Fit isoscape model with fixed effects of elevation and latitude
NqnFit <- isofit(data = rawWateragg, mean_model_fix = list(elev = TRUE, lat_abs = TRUE))
plot(NqnFit)
NqnFit
AIC(NqnFit$mean_fit)

# Load DEM raster
Elev <- rast("DEMNWPat180.tif")
plot(Elev)
extent(Elev)

# Prepare DEM raster for isoscape modeling
ElevNqn <- prepraster(raster = Elev, isofit = NqnFit, aggregation_factor = 4)
levelplot(ElevNqn,
          margin = FALSE,
          col.regions = viridisLite::turbo,
          main = "Structural raster")
ElevNqn
str(ElevNqn)

# Save processed DEM raster
writeRaster(ElevNqn,
            filename = "NNqn_DEMNew.tif")

# Generate isoscape from fitted model
NqnIsoscape <- isoscape(raster = ElevNqn, isofit = NqnFit)
NqnIsoscape$isoscapes
str(NqnIsoscape)

# Plot isoscape mean surface
plot(NqnIsoscape, y_title = list(which = TRUE, title = bquote(delta^18 * O)),
     sources = list(draw = FALSE), calibs = list(draw = FALSE), assigns = list(draw = FALSE))

# Save isoscape mean raster
writeRaster(NqnIsoscape$isoscapes$mean,
            filename = "NNqn_Isoscape3New.tif")

# Read human isotope datasets
humans <- read.csv("O18_humansC-HL.csv", h=TRUE)
humans
humans <- read.csv("O18_humansC-AqcoA.csv", h=TRUE)
humans
humans <- read.csv("O18_humansC-AqcoB.csv", h=TRUE)
humans
humans <- read.csv("O18_humansC-AqcoC.csv", h=TRUE)
humans

# Assign humans to isoscape
AssignedHumans <- isofind(data = humans,
                          isoscape = NqnIsoscape)
AssignedHumans

# Extract probability values for each human at their known location
extract(AssignedHumans$sample$pv[[1]], cbind(humans$long[1], humans$lat[1]))
Pvalues <- sapply(1:nrow(humans),
                  function(i) extract(AssignedHumans$sample$pv[[i]],
                                      cbind(humans$long[i], humans$lat[i])))
Pvalues
write.csv(Pvalues, "PvaluesAqcoC_2024New.csv")

# Identify samples with low probability (p <= 0.05)
humans$sample_ID[Pvalues <= 0.05]

# Plot assignment maps for selected individuals
plot(AssignedHumans, who = 1:6, sources = list(draw = FALSE),
     calibs = list(draw = FALSE),
     assigns = list(draw = TRUE, cex = 1.1, pch = 16, lwd = 1, col = "black"))
str(AssignedHumans)
AssignedHumans$sample$pv
humans
Pvalues
plot(AssignedHumans, who = 38:43, sources = list(draw = FALSE), calibs = list(draw = FALSE), assigns = list(draw = FALSE))
plot(AssignedHumans, who = 4:7, sources = list(draw = FALSE), calibs = list(draw = FALSE), assigns = list(draw = TRUE))
plot(AssignedHumans, who = 32:39, sources = list(draw = FALSE), calibs = list(draw = FALSE), assigns = list(draw = FALSE))
plot(AssignedHumans, who = 24:29, sources = list(draw = FALSE), calibs = list(draw = FALSE), assigns = list(draw = FALSE))
plot(AssignedHumans, who = "Aqcoind38")
plot(AssignedHumans, who = "Aqcoind31")
plot(AssignedHumans, who = "MILLAIN")
plot(AssignedHumans, who = "group")

# Inspect assignment object for specific sample
AssignedHumans$MILLAIN
str(AssignedHumans)

# Save isoscape and one individual probability raster
writeRaster(NqnIsoscape$isoscapes$mean,
            filename = "NqnIsoscape.tif",
            format = "GTiff",
            overwrite = TRUE,
            NAflag = -9999)
writeRaster(AssignedHumans$sample$pv$MILLAIN,
            filename = "MillainIsoscape.tif",
            format = "GTiff",
            overwrite = TRUE,
            NAflag = -9999)
