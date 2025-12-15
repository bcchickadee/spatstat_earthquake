# ============================================================================
# Script for the final project of the course "Statistical Analysis for Spatial Data"
# December 14, 2025. Sumin Yun, Seoul National University.
# ============================================================================

# Loading required libraries
library(tidyverse); 
library(sf); library(spatstat); library(terra); library(tidyterra); library(rnaturalearth)
library(magrittr); library(xtable)

# ============================================================================

# Importing earthquake dataset
df <- read_tsv("data/katalog_gempa_v2.tsv") %>%
  filter(magnitude >= 4.5, !is.na(latitude), !is.na(longitude))

# Outputting Table 1
df %>% 
  head() %>% 
  dplyr::select(datetime, latitude, longitude, magnitude, depth, location) %>% 
  mutate(datetime = datetime %>% strftime(format = "%Y-%m-%d %H:%M:%S")) %>% 
  xtable(caption = "Sample of earthquake data (only select columns are displayed)") %>% 
  print(comment = F)

# Transforming earthquake dataset to `sf` object
quake_sf <- st_as_sf(df, coords = c("longitude", "latitude"), crs = 4326) %>%
  st_transform(crs = 23845)

# Visualizing Figure 2
indonesia <- ne_countries(country = "Indonesia")
ggplot() +
  geom_sf(data = ne_countries(), fill = "grey") +
  geom_sf(data = indonesia, fill = "grey50") +
  geom_sf(data = quake_sf, mapping = aes(color = depth, size = magnitude), alpha = 0.1) +
  scale_color_viridis_c() +
  theme_bw() +
  labs(title = "Distribution of Earthquakes in Indonesia",
       x = "Longitude", y = "Latitude", size = "Magnitude", color = "Depth (km)") +
  coord_sf(xlim = c(95, 140), ylim = c(-12, 7))

# ============================================================================

# Importing fault dataset and transforming to `sf` object
fault_url <- "https://github.com/GEMScienceTools/gem-global-active-faults/raw/refs/heads/master/geopackage/gem_active_faults.gpkg"
temp_fault_gpkg <- tempfile(fileext = ".gpkg")
download.file(fault_url, destfile = temp_fault_gpkg, mode = "wb")
fault <- st_read(temp_fault_gpkg) %>% st_transform(crs = 4755) %>%
  filter(catalog_name == "EOS_SE_Asia") %>% 
  st_crop(xmin = min(df$longitude) - 3, xmax = max(df$longitude) + 3,
          ymin = min(df$latitude) - 3, ymax = max(df$latitude) + 3) %>% 
  mutate(average_dip = average_dip %>% parse_number(),
         average_rake = average_rake %>% parse_number())

fault$geom <- fault$geom %>% st_transform(crs = 23845)
file.remove(temp_fault_gpkg)

# Visualizing Figure 3
ggplot() +
  geom_sf(data = ne_countries(), fill = "grey") +
  geom_sf(data = indonesia, fill = "grey50") +
  geom_sf(data = fault, mapping = aes(color = average_dip), linewidth = 1.5) +
  scale_color_viridis_c() +
  theme_bw() +
  labs(title = "Location of Faults in Indonesia",
       x = "Longitude", y = "Latitude", color = "Average Dip (°)") +
  coord_sf(xlim = c(92, 142), ylim = c(-12, 8))

# ============================================================================

# Importing subduction zone dataset
subduction_files <- c(
  (c("sum", "sul", "hal") %>%
     paste0("data/slab2/", ., "_slab2_dep_02.23.18.grd")),
  (c("phi", "png") %>%
     paste0("data/slab2/", ., "_slab2_dep_02.26.18.grd"))
  )
subduction <- sprc(lapply(subduction_files, rast)) %>%
  mosaic(fun = "max") %>% project("EPSG:23845")

# Visualizing Figure 4
ggplot() +
  geom_sf(data = ne_countries(), fill = "grey") +
  geom_sf(data = indonesia, fill = "grey50") +
  geom_spatraster(data = subduction) +
  scale_fill_viridis_c(option = "plasma", na.value = NA) +
  theme_bw() +
  labs(title = "Depth Distribution of Subduction Zones near Indonesia",
       x = "Longitude", y = "Latitude", fill = "Depth (km)") +
  coord_sf(xlim = c(85, 152), ylim = c(-17, 30))

# ============================================================================

# Visualizing Figure 5

fault_dist <- fault$geom %>% as.psp() %>% distfun() %>%
  as.im() %>% as.data.frame() %>%  rast(crs = "EPSG:23845")
ggplot() +
  geom_spatraster(data = fault_dist) +
  scale_fill_viridis_c(option = "plasma") +
  labs(title = "Distance to nearest active fault in Indonesia",
       x = "Longitude", y = "Latitude", fill = "Distance (m)") +
  theme_bw()

# Visualizing Figure 6

fault_dist %>%
  hist(breaks = 100,
       main = "Histogram of distance to nearest fault in Indonesia",
       xlab = "Distance (m)")

# Visualizing Figure 7

subduction %>% 
  hist(breaks = 100,
       main = "Histogram of depth of subduction zone in Indonesia",
       xlab = "Distance (m)")

# ============================================================================

# Making `ppp` object from the earthquake `sf` object
quake_coords <- st_coordinates(quake_sf)
quake_window <- owin(
  c(quake_coords[,1] %>% min(), quake_coords[,1] %>% max()),
  c(quake_coords[,2] %>% min(), quake_coords[,2] %>% max())
)

quake_ppp <- ppp(x = quake_coords[,1],
                 y = quake_coords[,2],
                 window = quake_window,
                 unitname = c("degree", "degrees"))

if(any(duplicated(quake_ppp))) {
  quake_ppp <- rjitter(quake_ppp, retry = T, nsim = 1, drop = T)
}

# Counting quadrats and visualizing Figure 8
quake_count <- quadratcount(quake_ppp, nx = 20, ny = 10)
quake_count %>% as.data.frame() %>% tibble() %>% 
  separate(col = "x", sep = ",", into = c("x_lower", "x_upper")) %>% 
  separate(col = "y", sep = ",", into = c("y_lower", "y_upper")) %>% 
  mutate(x_lower = parse_number(x_lower),
         x_upper = parse_number(x_upper),
         y_lower = parse_number(y_lower),
         y_upper = parse_number(y_upper),
         x = (x_lower + x_upper) / 2,
         y = (y_lower + y_upper) / 2) %>% 
  select(x, y, Freq) %>% 
  ggplot(mapping = aes(x = x, y = y, fill = Freq)) +
  geom_tile() +
  scale_fill_viridis_c(option = "plasma", trans = "sqrt") +
  geom_text(mapping = aes(label = Freq), angle = 90) +
  theme_bw() +
  labs(title = "Quadrat Count of Earthquakes in Indonesia", fill = "Count") +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank())

# Testing CSR with `quadrat.test`
quadrat.test(quake_ppp, nx = 20, ny = 10)

# ============================================================================

# Kernel density estimation
sigma_cv <- bw.diggle(quake_ppp)
lambda_cv <- density(quake_ppp, sigma = sigma_cv, positive = T,
                     kernel = "gaussian", diggle = T)

# Visualizing Figure 9
ggplot() +
  geom_spatraster(data = rast(as.data.frame(lambda_cv),
                              crs = "EPSG:23845")) +
  scale_fill_viridis_c(option = "plasma", na.value = NA, trans = "sqrt") +
  theme_bw() +
  labs(title = "Estimated Intensity of Earthquakes in Indonesia",
       x = "Longitude", y = "Latitude", fill = "Intensity")

# ============================================================================

# Getting Ripley's K and L-functions and visualizing Figure 10
par(mfrow = c(1,2))
Kest(quake_ppp) %>% plot(main = "Estimated K-function", xlab = "Distance (m)")
Lest(quake_ppp) %>% plot(main = "Estimated L-function", xlab = "Distance (m)")

# ============================================================================

# Inhomogeneous K-function
K_in <- Kinhom(quake_ppp, lambda = lambda_cv, correction = "translate")
env_K <- envelope(quake_ppp, Kinhom, lambda = lambda_cv,
                  nsim = 30, global = TRUE)


# Visualizing Figure 11
plot(env_K, sqrt(./pi) ~ r, main="Inhomogeneous L-Function with Envelopes",
     xlab = "Distance (m)")

# ============================================================================

# Fitting LGCP Model
fault_psp <- fault$geom %>% 
  as.psp()

D_faults <- distfun(fault_psp)

subduction_df <- subduction %>% as.data.frame(xy = T)
colnames(subduction_df) <- c("x", "y", "depth")
subduction_im <- as.im(subduction_df)

fit_lgcp <- kppm(
  quake_ppp ~ D_faults + subduction_depth,
  clusters = "LGCP",
  model = "exponential",
  covariates = list(D_faults = D_faults, subduction_depth = subduction_im)
)

# Outputting Table 2
lgcp_coef <- fit_lgcp %>% summary() %>% coef()

lgcp_coef %>%
  xtable(caption = "Parameters of fitted coefficients of LGCP model",
         digits = c(-4, -4, -4, -4, -4, 4, 4)) %>% 
  print(comment = F)

# Visualizing Figure 12
lambda_lgcp <- predict(fit_lgcp)
ggplot() +
  geom_sf(data = ne_countries(), fill = "grey") +
  geom_sf(data = indonesia, fill = "grey50") +
  geom_spatraster(data = rast(as.data.frame(lambda_lgcp), crs = "EPSG:23845")) +
  scale_fill_viridis_c(option = "plasma", na.value = NA) +
  theme_bw() +
  labs(title = "Predicted Intensity from LGCP Model",
       x = "Longitude", y = "Latitude", fill = "Intensity") +
  coord_sf(xlim = c(92, 145), ylim = c(-14, 10))

# Visualizing Figure 13
par(mfrow = c(1, 2))
plot(fit_lgcp, what="statistic", pause=FALSE,
     main = "Fitted Inhomogeneous\nK-Function using LGCP",
     xlab = "Distance (m)")
plot(fit_lgcp, sqrt(./pi) ~ r, what="statistic", pause=FALSE,
     main = "Fitted Inhomogeneous\nL-Function using LGCP",
     xlab = "Distance (m)")

# Visualizing Figure 14
plot(fit_lgcp,
     sqrt(./pi) - r ~ r,
     what="statistic", pause=FALSE,
     main = "Fitted Inhomogeneous L-Function using LGCP",
     xlab = "Distance (m)",
     ylab = expression(sqrt(K[inhom](r)/pi) - r),
     legend = F)

legend(
  "topright",
  inset = 0.02,
  legend = c(
    expression(sqrt(K[fit](r)/pi) - r),
    expression(sqrt(K[inhom]^iso(r)/pi) - r),
    expression(sqrt(K[inhom]^pois(r)/pi) - r)
  ),
  col = c("black", "red", "darkgreen"),
  lty = c(1, 2, 3),
  lwd = c(2, 2, 1),
  cex = 0.9,
  bty = "o"
)
