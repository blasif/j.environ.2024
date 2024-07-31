# New DATASET application

# rm(list = ls())
source('00_system.R')

ifelse(length(which(installed.packages()[,1] == 'ncdf4')) == 0 , 
       install.packages('ncdf4'), print('ncdf4 installed'))

path <- paste0('Datasets/')

file_prec <- 'precipitation_monthly-mean_era5-to-1km_1979-2018-mean_v1.0.nc'
file_wind <- 'wind-speed_monthly-mean_era5-to-1km_1979-2018-mean_v1.0.nc'
file_cc <- 'cloud-cover_monthly-mean_era5-to-1km_1979-2018-mean_v1.0.nc'
fil_BIO_04 <- 'BIO04_era5-to-1km_1979-2018-mean_v1.0.nc'
file_BIO_15 <- 'BIO15_era5-to-1km_1979-2018-mean_v1.0.nc'
file_wind_meridional <- 'meridional-wind-speed_monthly-mean_era5-to-1km_1979-2018-median_v1.0.nc'

nc_prec <- nc_open(paste0(path, file_prec), auto_GMT = T)
lati <- ncvar_get(nc_prec, 'latitude')
long <- ncvar_get(nc_prec, 'longitude')
prec <- ncvar_get(nc_prec, 'precipitation_monthly-mean')

full_df <- cbind.data.frame(expand.grid('long' = long, 'lati' = lati),
                            'prec' = c(prec[, , 1]))

rm(prec, long, lati)

nc_wind <- nc_open(paste0(path, file_wind), auto_GMT = T)
nc_cc <- nc_open(paste0(path, file_cc), auto_GMT = T)
nc_BIO_04 <- nc_open(paste0(path, fil_BIO_04), auto_GMT = T)
nc_BIO_15 <- nc_open(paste0(path, file_BIO_15), auto_GMT = T)
nc_wind_meri <- nc_open(paste0(path, file_wind_meridional), auto_GMT = T)

wind <- ncvar_get(nc_wind, 'wind-speed')

full_df <- cbind.data.frame(full_df,
                            'wind' = c(wind[, , 1]))

rm(wind)

cloud_cover <- ncvar_get(nc_cc, 'cloud-cover')

full_df <- cbind.data.frame(full_df,
                            'cloud_c' = c(cloud_cover[,,1]))

rm(cloud_cover)

BIO_04 <- ncvar_get(nc_BIO_04, 'BIO04')

full_df <- cbind.data.frame(full_df,
                            'BIO04' = c(BIO_04)
)

rm(BIO_04)

BIO_15 <- ncvar_get(nc_BIO_15, 'BIO15')

full_df <- cbind.data.frame(full_df,
                            'BIO15' = c(BIO_15)
)

rm(BIO_15)

wind_mer <- ncvar_get(nc_wind_meri, 'meridional-wind-speed')

full_df <- cbind.data.frame(full_df,
                            'new_wind' = c(wind_mer[,,1])
)

rm(wind_mer)

full_df$index <- 1:dim(full_df)[1]

# first filtering
index_to_swiss <- which(full_df$long < 10.75 & full_df$long > 5 & 
                          full_df$lati < 48 & full_df$lati > 45.75)

swiss_df <- full_df[index_to_swiss, ]

poly_sf <- st_as_sf(swiss_df, coords = c('long','lati'), crs = st_crs("+proj=longlat +datum=WGS84 +no_defs"))

country_boundary <- st_as_sf(geodata::gadm("Switzerland", level=1, path = "."))

#

country_boundary_transf <- st_transform(country_boundary, crs = st_crs("+proj=longlat +datum=WGS84 +no_defs"))

intersection_step <- st_intersection(poly_sf, country_boundary_transf)

df_step <- cbind.data.frame(as.data.frame(intersection_step), 
                            st_coordinates(intersection_step))
df_step <- df_step[, -c(8:19)]
colnames(df_step)[8:9] <- c("long", "lati")

# TO get the elevation

swiss_elev <- df_step
coordinates(swiss_elev) <- c('long', 'lati')
proj4string(swiss_elev) <- ("+proj=longlat +datum=WGS84 +no_defs")

eee <- st_as_sf(swiss_elev)

aaa <- elevatr::get_elev_point(eee, units = 'meters', 
                               src = 'aws', z = 7)

swiss_finalll <-as.data.frame(aaa)

df_step$elevation <- swiss_finalll$elevation

df_step <- df_step[,c(8,9,1,2,3,4,5,6,10)]

df_step$log_elevation <- log(df_step$elevation)
df_step$log_wind <- log(df_step$wind)
df_step$log_cloud <- log(df_step$cloud_c)
df_step$log_lati <- log(df_step$lati)
df_step$log_long <- log(df_step$long)
df_step$prec <- df_step$prec * 60 * 60 * 24 * 1000
# Create subsets

set.seed(181222)

# Create clusters

if(T){
  
  # Cluster based on longitude, latitude, and cloud coverage
  
  to_test <-as.data.frame(scale(df_step[,c(1,2,5)]))
  
  a_ver_che <- kmeans(to_test,centers = 475,iter.max = 50)
  
  to_check <- sample(500,50)
  
  striped_data_b <- which(a_ver_che$cluster %in% to_check)

}

# Create stripes

if(T){
  
  angle_deg <- 10
  
  # Function to rotate coordinates
  rotate_coords <- function(lat, long, angle_deg) {
    angle_rad <- angle_deg * pi / 180
    lat_rot <- lat * cos(angle_rad) - long * sin(angle_rad)
    long_rot <- lat * sin(angle_rad) + long * cos(angle_rad)
    return(list(lat_rot = lat_rot, long_rot = long_rot))
  }
  
  # Rotate the DataFrame coordinates
  coords_rotated <- rotate_coords(df_step$lati, df_step$long, angle_deg)
  df_step$lat_rot <- coords_rotated$lat_rot
  df_step$long_rot <- coords_rotated$long_rot
  
  strip_width <- 0.005  # width of each stripe and gap
  strips <- cut(df_step$long_rot, breaks=seq(min(df_step$long_rot), max(df_step$long_rot), by=strip_width))
  
  strip_width <- 0.01  # width of each stripe and gap
  strips_2 <- cut(df_step$long_rot, breaks=seq(min(df_step$long_rot), max(df_step$long_rot), by=strip_width))
  
  strip_width <- 0.05  # width of each stripe and gap
  strips_3 <- cut(df_step$long_rot, breaks=seq(min(df_step$long_rot), max(df_step$long_rot), by=strip_width))
  
  strip_width <- 0.1  # width of each stripe and gap
  strips_4 <- cut(df_step$long_rot, breaks=seq(min(df_step$long_rot), max(df_step$long_rot), by=strip_width))
  
  striped_data <- which(strips %in% levels(strips)[sample(length(levels(strips)),40)] |
                          strips_2 %in% levels(strips_2)[sample(length(levels(strips_2)),20)] |
                          strips_3 %in% levels(strips_3)[sample(length(levels(strips_3)),2)] |
                          strips_4 %in% levels(strips_4)[sample(length(levels(strips_4)),2)])
  
  striped_data_f <- unique(c(striped_data, striped_data_b))
  
  print(length(striped_data_f))
  
  aaa <- sample(69965, 20000 - length(striped_data_f))
  
  striped_data_ff <- unique(c(striped_data_f, aaa))
  
  all_dfs <- list()
  
  df_step <- df_step[, -c(15, 16)]
  
  all_dfs[[1]] <- df_step[striped_data_ff, ]
  
  # training_sample
  
  df_remaining_two <- df_step[-striped_data_ff, ]
  
  # random_sample <- sample(dim(df_remaining_two)[1], size = 500)
  
  seq_sample <- floor(seq(1, dim(df_remaining_two)[1], length.out = 500))
  all_dfs[[2]] <- df_remaining_two[seq_sample, ]
  
  df_reamining_three <- df_remaining_two[-seq_sample, ]
  
  # for tapering
  
  seq_sample_sparse <- floor(seq(1,dim(df_reamining_three)[1], length.out = 10000))
  all_dfs[[3]] <- df_reamining_three[seq_sample_sparse, ]
  
  mean(all_dfs[[1]]$prec) - mean(all_dfs[[2]]$prec)
  
  all_dfs[[4]] <- df_step
  
  save(all_dfs, file = 'Datasets/app_dataset.RData')
  
}
