threed_data <- read.csv("Data/Old data/Jacob/3D_MECP2_6WO_combined_dataset.csv")
threed_data <- threed_data %>%
  select(Filename, Mean, CX..pix., CY..pix., CZ..pix.) %>%
  rename_all(tolower) %>%
  rename(
    x = cx..pix.,
    y = cy..pix.,
    z = cz..pix.
  ) %>%
  mutate(hemisphere = str_extract(filename, "(LH|RH)")) %>%
  mutate(condition = str_extract(filename, "(NW|NH|SW|SH)")) %>%
  mutate(mecp2_p = "P") %>%
  mutate(id = 1:nrow(.))

threed_data <- threed_data %>% group_by(filename)
all_imgs <- threed_data %>% group_split(threed_data)

# Function to calculate Euclidean distance between two 3D points
euclidean_distance_3d <- function(x1, y1, z1, x2, y2, z2) {
  sqrt((x2 - x1)^2 + (y2 - y1)^2 + (z2 - z1)^2)
}

# Calculate nearest neighbor for each row
for (i in 1:nrow(df)) {
  observation <- df[i, ]
  other_points <- df[-i, ]
  distances <- euclidean_distance_3d(
    observation$x, observation$y, observation$z,
    other_points$x, other_points$y, other_points$z
  )
  nearest_index <- which.min(distances)
  nearest_point <- other_points[nearest_index, ]
  
  print(paste("Nearest neighbor for row", i, "is row", nearest_index + 1))
  print(nearest_point)
}


# Calculate nearest neighbor for each row on actual data
for (i in 1:length(all_imgs)) {
  df <- all_imgs[[i]]
  for ( r in 1:nrow(df)){
    observation <- df[r, ]
    other_points <- df[-r, ]
    distances <- euclidean_distance_3d(
      observation$x, observation$y, observation$z,
      other_points$x, other_points$y, other_points$z
    )
    nearest_index <- which.min(distances)
    nearest_point <- other_points[nearest_index, ]
    
    print(paste("Nearest neighbor for row", r, "is row", nearest_index))
    print(nearest_point)
  }
}


# For multiple neighbors
k_neighbors <- 5
for (i in 1:length(all_imgs)) {
  for (r in 1:nrow(df)) {
    observation <- df[r, ]
    other_points <- df[-r, ]
    distances <- euclidean_distance_3d(
      observation$x, observation$y, observation$z,
      other_points$x, other_points$y, other_points$z
    )
    k_nearest_indices <- order(distances)[1:k_neighbors]
    k_nearest_points <- other_points[k_nearest_indices, ]

    print(paste("Nearest", k_neighbors, "neighbors for row", r, ":"))
    print(k_nearest_points)
  }
}

