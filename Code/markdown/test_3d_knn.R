threed_data <- read.csv("Data/Old data/Jacob/3D_MECP2_6WO_combined_dataset.csv")
threed_data <- threed_data %>%
  select(Filename, Mean, CX..pix., CY..pix., CZ..pix., MECP2) %>%
  rename(mecp2_p = MECP2) %>%
  rename_all(tolower) %>%
  rename(
    x = cx..pix.,
    y = cy..pix.,
    z = cz..pix.,
  ) %>%
  mutate(hemisphere = str_extract(filename, "(LH|RH)")) %>%
  mutate(condition = str_extract(filename, "(NW|NH|SW|SH)")) %>%
  filter(mecp2_p == "P") %>%
  mutate(id = 1:nrow(.))

threed_data <- threed_data %>% group_by(filename)
all_imgs <- threed_data %>% group_split(threed_data)



# Function to calculate Euclidean distance between two 3D points
euclidean_distance_3d <- function(x1, y1, z1, x2, y2, z2) {
  sqrt((x2 - x1)^2 + (y2 - y1)^2 + (z2 - z1)^2)
}

all_data <- list()
all_p_scores <- list()
counter <- 1
apply_intensity <- TRUE
p_only <- TRUE
calc_mean_intensity <- TRUE
neighbor_means <- c()
# For multiple neighbors
k_neighbors <- 5
for (i in 1:length(all_imgs)) {
  df <- all_imgs[[i]]
  df$knn_5_label <- rep(0, nrow(df))
  for (r in 1:nrow(df)) {
    observation <- df[r, ]
    other_points <- df[-r, ]
    distances <- euclidean_distance_3d(
      observation$x, observation$y, observation$z,
      other_points$x, other_points$y, other_points$z
    )
    
    
    k_nearest_indices <- order(distances)[1:k_neighbors]
    k_nearest_points <- other_points[k_nearest_indices, ]
    
    intensity_info <- other_points[k_nearest_indices, ]
    add_p_score <- sum((intensity_info$mecp2_p == "P")) / k_neighbors
    observation$knn_5_label <- add_p_score
    all_p_scores[[counter]] <- add_p_score
    
    # Calculate distances with intensity information
    if (apply_intensity) {
      distances <- distances[k_nearest_indices]
      distances_updated <- rescale(distances * intensity_info$mean, to = c(0, 1))
      distances_updated <- sort(distances_updated, decreasing = FALSE)

      # Calculate updated distances with intensity information for just P neighbors
      if (p_only) {
        intensity_info <- filter(intensity_info, mecp2_p == "P")
        if (nrow(intensity_info) == 0) {
          writeLines("This sample has only N neighbors.\nIt is being excluded from this analysis.")
        }
      }
      if (calc_mean_intensity) {
        # intensity_info$positive_neighborhood_mean <- rep(0, nrow(intensity_info))
        df$positive_neighborhood_mean <- rep(0, nrow(df))
        mean_of_neighbors <- mean(intensity_info$mean)
        neighbor_means <- c(neighbor_means, mean_of_neighbors)
        observation$positive_neighborhood_mean <- mean_of_neighbors
        all_data[[counter]] <- observation
        counter <- counter + 1
      }
    }

    # print(paste("Nearest", k_neighbors, "neighbors for row", r, ":"))
    # print(k_nearest_points)
  }
}

test_df <- bind_rows(all_data)

# + Neighborhood mean
# Now trying the neighborhood plots again
all_imgs <- list()
counter <- 1
mean_neighbor_df_simplified <- test_df
mean_neighbor_df_simplified <- mean_neighbor_df_simplified %>%
  mutate(filename = match(filename, unique(filename)))

for (i in unique(mean_neighbor_df_simplified$condition)) {
  current_cond <- filter(mean_neighbor_df_simplified, condition == i)
  current_cond$filename <- as.factor(current_cond$filename)
  
  # Test for normality using Shapiro-Wilk test
  shapiro_test <- shapiro.test(current_cond$positive_neighborhood_mean)
  if (shapiro_test$p.value >= 0.05) {
    # Data is normally distributed
    print("Data is normally distributed")
    
    # Perform Levene's test for homogeneity of variances
    levene_test <- leveneTest(positive_neighborhood_mean ~ filename, data = current_cond)
    print(levene_test)
    
    p_value <- levene_test$"Pr(>F)"[1] # Get the p-value from the test
    p_value_text <- ifelse(p_value < 0.05, paste0("p = ", format(p_value, scientific = FALSE, digits = 3)), "NS")
  } else {
    # Data is not normally distributed
    print("Data is not normally distributed")
    
    # Perform Fligner-Killeen test for homogeneity of variances
    fligner_test <- fligner.test(positive_neighborhood_mean ~ filename, data = current_cond)
    print(fligner_test)
    
    
    kw_test <- kruskal.test(positive_neighborhood_mean ~ filename, data = current_cond)
    print(kw_test)
    
    p_value <- fligner_test$p.value
    p_value_kw <- kw_test$p.value
    p_value_text <- ifelse(p_value < 0.05, paste0("p = ", format(p_value, scientific = TRUE, digits = 3)), "NS")
    p_value_text <- ifelse(p_value_kw < 0.05, paste0("p = ", format(p_value_kw, scientific = TRUE, digits = 3)), "NS")
  }
  
  
  
  p1 <- ggplot(aes(x = filename, y = positive_neighborhood_mean), data = current_cond) +
    geom_violin(draw_quantiles = c(0.25, 0.50, 0.75)) +
    geom_sina(aes(colour = filename))+
    theme(
      panel.background = element_blank(),
      axis.title = element_text(size = 18),
      axis.text = element_text(size = 16),
      axis.ticks.length = unit(0.2, units = "cm"),
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
      axis.line = element_line(color = "black", lineend = "round", linewidth = 1.25),
      legend.title = element_blank(),
      legend.text = element_blank(),
      legend.key = element_blank(),
      legend.position = "none"
    ) +
    ggtitle(paste0("K = 5 | ", current_cond$condition, " | ", p_value_text)) +
    xlab("Image") +
    ylab("+ Neighborhood Mean") +
    coord_flip()
  all_imgs[[counter]] <- p1
  counter <- counter + 1
}



p2 <- ggarrange(all_imgs[[4]],
                all_imgs[[3]],
                all_imgs[[2]],
                all_imgs[[1]],
                ncol = 2,
                nrow = 2,
                align = "hv",
                labels = "AUTO"
)





# Mean Intensity of Cells
# Now trying the neighborhood plots again
all_imgs <- list()
counter <- 1
mean_neighbor_df_simplified <- test_df
mean_neighbor_df_simplified <- mean_neighbor_df_simplified %>%
  mutate(filename = match(filename, unique(filename)))

for (i in unique(mean_neighbor_df_simplified$condition)) {
  current_cond <- filter(mean_neighbor_df_simplified, condition == i)
  current_cond$filename <- as.factor(current_cond$filename)
  
  # Test for normality using Shapiro-Wilk test
  shapiro_test <- shapiro.test(current_cond$mean)
  if (shapiro_test$p.value >= 0.05) {
    # Data is normally distributed
    print("Data is normally distributed")
    
    # Perform Levene's test for homogeneity of variances
    levene_test <- leveneTest(mean ~ filename, data = current_cond)
    print(levene_test)
    
    p_value <- levene_test$"Pr(>F)"[1] # Get the p-value from the test
    p_value_text <- ifelse(p_value < 0.05, paste0("p = ", format(p_value, scientific = FALSE, digits = 3)), "NS")
  } else {
    # Data is not normally distributed
    print("Data is not normally distributed")
    
    # Perform Fligner-Killeen test for homogeneity of variances
    fligner_test <- fligner.test(mean ~ filename, data = current_cond)
    print(fligner_test)
    
    
    kw_test <- kruskal.test(mean ~ filename, data = current_cond)
    print(kw_test)
    
    p_value <- fligner_test$p.value
    p_value_kw <- kw_test$p.value
    p_value_text <- ifelse(p_value < 0.05, paste0("p = ", format(p_value, scientific = TRUE, digits = 3)), "NS")
    # p_value_text <- ifelse(p_value_kw < 0.05, paste0("p = ", format(p_value_kw, scientific = TRUE, digits = 3)), "NS")
  }
  
  
  
  p1 <- ggplot(aes(x = filename, y = mean), data = current_cond) +
    geom_violin(draw_quantiles = c(0.25, 0.50, 0.75)) +
    geom_sina(aes(colour = filename))+
    theme(
      panel.background = element_blank(),
      axis.title = element_text(size = 18),
      axis.text = element_text(size = 16),
      axis.ticks.length = unit(0.2, units = "cm"),
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
      axis.line = element_line(color = "black", lineend = "round", linewidth = 1.25),
      legend.title = element_blank(),
      legend.text = element_blank(),
      legend.key = element_blank(),
      legend.position = "none"
    ) +
    ggtitle(paste0("K = 5 | ", current_cond$condition, " | ", p_value_text)) +
    xlab("Image") +
    ylab("Mean Intensity of Cells") +
    coord_flip()
  all_imgs[[counter]] <- p1
  counter <- counter + 1
}



p2 <- ggarrange(all_imgs[[4]],
                all_imgs[[3]],
                all_imgs[[2]],
                all_imgs[[1]],
                ncol = 2,
                nrow = 2,
                align = "hv",
                labels = "AUTO"
)


# Now doing scatter plots
plots <- list()
counter <- 1

for (c in unique(mean_neighbor_df_simplified$condition)) {
  current_cond <- filter(mean_neighbor_df_simplified, condition == c)
  current_cond$filename <- as.factor(current_cond$filename)
  p1 <- ggscatter(
    data = current_cond,
    x = "mean", y = "positive_neighborhood_mean",
    conf.int = TRUE,
    cor.coef = TRUE,
    cor.coeff.args = list(
      method = "spearman",
      label.x = 1000, label.y = 3000, label.sep = ",",
      cor.coef.name = "rho"
    ),
    parse = TRUE,
    cor.coef.size = 10,
    shape = 21,
    size = 2.5,
    add = fit_line,
    add.params = list(
      color = "black",
      fill = "grey"
    ),
    repel = TRUE,
    label.rectangle = TRUE,
    color = "filename",
    fill = "filename",
    alpha = 1,
    cor.method = stat_used,
    xlab = "Mean Intensity of Cells",
    ylab = "Mean MECP2+\nNeighborhood Intensity",
    ellipse = FALSE,
    mean.point = FALSE
  )
  
  # p2 for the single annotation context
  p2 <- ggpar(p1,
              title = paste0("K = ", k_n, " | ", c, " | 3D "),
              font.main = c(20, "bold"),
              font.x = c(18, "bold"),
              font.y = c(18, "bold"),
              ggtheme = theme(plot.title = element_text(hjust = 0.5))
  )
  
  p2 <- p2 +
    theme(
      axis.line = element_line(linewidth = 1.25, lineend = "round"),
      axis.text = element_text(size = 18),
      axis.ticks = element_line(linewidth = 0.2),
      plot.subtitle = element_text(face = "italic"),
      legend.title = element_text(size = 16, face = "bold"),
      legend.text = element_text(size = 14),
      legend.position = "bottom"
    )+
    labs(color = "Image", fill = "Image")
  
  # Change to p2 for single annotation case, p5 otherwise
  plots[[paste0(c)]] <- p2
  counter <- counter + 1
}

# Combo plot
combo_plot <- ggarrange(
  plots$NW,
  plots$NH,
  plots$SW,
  plots$SH,
  ncol = 2,
  nrow = 2,
  labels = "AUTO",
  align = "hv",
  legend = "right",
  common.legend = FALSE,
  font.label = list(
    size = 18,
    face = "bold"
  )
)


# 12 week dataset
twelve_df <- read.csv("Data/Old data/Jacob/12WOcombined_output_MECP2.csv")
twelve_df <- twelve_df %>%
  select(Image, CX..pix., CY..pix., CZ..pix., Condition, Hemishphere, MECP2, Mean) %>%
  rename(x = CX..pix.,
         y = CY..pix.,
         z = CZ..pix.,
         mecp2_p = MECP2) %>%
  rename_all(tolower) %>%
  filter(mecp2_p == "P")

twelve_df <- twelve_df %>% group_by(image)
all_imgs <- twelve_df %>% group_split(twelve_df)


all_data <- list()
all_p_scores <- list()
counter <- 1
apply_intensity <- TRUE
p_only <- TRUE
calc_mean_intensity <- TRUE
neighbor_means <- c()

# For multiple neighbors
k_neighbors <- 5
for (i in 1:length(all_imgs)) {
  df <- all_imgs[[i]]
  df$knn_5_label <- rep(0, nrow(df))
  for (r in 1:nrow(df)) {
    observation <- df[r, ]
    other_points <- df[-r, ]
    distances <- euclidean_distance_3d(
      observation$x, observation$y, observation$z,
      other_points$x, other_points$y, other_points$z
    )
    
    
    k_nearest_indices <- order(distances)[1:k_neighbors]
    k_nearest_points <- other_points[k_nearest_indices, ]
    
    intensity_info <- other_points[k_nearest_indices, ]
    add_p_score <- sum((intensity_info$mecp2_p == "P")) / k_neighbors
    observation$knn_5_label <- add_p_score
    all_p_scores[[counter]] <- add_p_score
    
    # Calculate distances with intensity information
    if (apply_intensity) {
      distances <- distances[k_nearest_indices]
      distances_updated <- rescale(distances * intensity_info$mean, to = c(0, 1))
      distances_updated <- sort(distances_updated, decreasing = FALSE)
      
      # Calculate updated distances with intensity information for just P neighbors
      if (p_only) {
        intensity_info <- filter(intensity_info, mecp2_p == "P")
        if (nrow(intensity_info) == 0) {
          writeLines("This sample has only N neighbors.\nIt is being excluded from this analysis.")
          next
        }
      }
      if (calc_mean_intensity) {
        # intensity_info$positive_neighborhood_mean <- rep(0, nrow(intensity_info))
        df$positive_neighborhood_mean <- rep(0, nrow(df))
        mean_of_neighbors <- mean(intensity_info$mean)
        neighbor_means <- c(neighbor_means, mean_of_neighbors)
        observation$positive_neighborhood_mean <- mean_of_neighbors
        all_data[[counter]] <- observation
        counter <- counter + 1
      }
    }
    
    # print(paste("Nearest", k_neighbors, "neighbors for row", r, ":"))
    # print(k_nearest_points)
  }
}

test_df <- bind_rows(all_data)
test_df <- test_df %>%
  rename(filename = image)

# + Neighborhood mean
# Now trying the neighborhood plots again
all_imgs <- list()
counter <- 1
mean_neighbor_df_simplified <- test_df
mean_neighbor_df_simplified <- mean_neighbor_df_simplified %>%
  mutate(filename = match(filename, unique(filename)))

for (i in unique(mean_neighbor_df_simplified$condition)) {
  current_cond <- filter(mean_neighbor_df_simplified, condition == i)
  current_cond$filename <- as.factor(current_cond$filename)
  
  # Test for normality using Shapiro-Wilk test
  shapiro_test <- shapiro.test(current_cond$positive_neighborhood_mean)
  if (shapiro_test$p.value >= 0.05) {
    # Data is normally distributed
    print("Data is normally distributed")
    
    # Perform Levene's test for homogeneity of variances
    levene_test <- leveneTest(positive_neighborhood_mean ~ filename, data = current_cond)
    print(levene_test)
    
    p_value <- levene_test$"Pr(>F)"[1] # Get the p-value from the test
    p_value_text <- ifelse(p_value < 0.05, paste0("p = ", format(p_value, scientific = FALSE, digits = 3)), "NS")
  } else {
    # Data is not normally distributed
    print("Data is not normally distributed")
    
    # Perform Fligner-Killeen test for homogeneity of variances
    fligner_test <- fligner.test(positive_neighborhood_mean ~ filename, data = current_cond)
    print(fligner_test)
    
    
    kw_test <- kruskal.test(positive_neighborhood_mean ~ filename, data = current_cond)
    print(kw_test)
    
    p_value <- fligner_test$p.value
    p_value_kw <- kw_test$p.value
    p_value_text <- ifelse(p_value < 0.05, paste0("p = ", format(p_value, scientific = TRUE, digits = 3)), "NS")
    p_value_text <- ifelse(p_value_kw < 0.05, paste0("p = ", format(p_value_kw, scientific = TRUE, digits = 3)), "NS")
  }
  
  
  
  p1 <- ggplot(aes(x = filename, y = positive_neighborhood_mean), data = current_cond) +
    geom_violin(draw_quantiles = c(0.25, 0.50, 0.75)) +
    geom_sina(aes(colour = filename))+
    theme(
      panel.background = element_blank(),
      axis.title = element_text(size = 18),
      axis.text = element_text(size = 16),
      axis.ticks.length = unit(0.2, units = "cm"),
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
      axis.line = element_line(color = "black", lineend = "round", linewidth = 1.25),
      legend.title = element_blank(),
      legend.text = element_blank(),
      legend.key = element_blank(),
      legend.position = "none"
    ) +
    ggtitle(paste0("K = 5 | ", current_cond$condition, " | ", p_value_text)) +
    xlab("Image") +
    ylab("+ Neighborhood Mean") +
    coord_flip()
  all_imgs[[counter]] <- p1
  counter <- counter + 1
}



p2 <- ggarrange(all_imgs[[4]],
                all_imgs[[3]],
                all_imgs[[2]],
                all_imgs[[1]],
                ncol = 2,
                nrow = 2,
                align = "hv",
                labels = "AUTO"
)





# Mean Intensity of Cells
# Now trying the neighborhood plots again
all_imgs <- list()
counter <- 1
mean_neighbor_df_simplified <- test_df
mean_neighbor_df_simplified <- mean_neighbor_df_simplified %>%
  mutate(filename = match(filename, unique(filename)))

for (i in unique(mean_neighbor_df_simplified$condition)) {
  current_cond <- filter(mean_neighbor_df_simplified, condition == i)
  current_cond$filename <- as.factor(current_cond$filename)
  
  # Test for normality using Shapiro-Wilk test
  shapiro_test <- shapiro.test(current_cond$mean)
  if (shapiro_test$p.value >= 0.05) {
    # Data is normally distributed
    print("Data is normally distributed")
    
    # Perform Levene's test for homogeneity of variances
    levene_test <- leveneTest(mean ~ filename, data = current_cond)
    print(levene_test)
    
    p_value <- levene_test$"Pr(>F)"[1] # Get the p-value from the test
    p_value_text <- ifelse(p_value < 0.05, paste0("p = ", format(p_value, scientific = FALSE, digits = 3)), "NS")
  } else {
    # Data is not normally distributed
    print("Data is not normally distributed")
    
    # Perform Fligner-Killeen test for homogeneity of variances
    fligner_test <- fligner.test(mean ~ filename, data = current_cond)
    print(fligner_test)
    
    
    kw_test <- kruskal.test(mean ~ filename, data = current_cond)
    print(kw_test)
    
    p_value <- fligner_test$p.value
    p_value_kw <- kw_test$p.value
    p_value_text <- ifelse(p_value < 0.05, paste0("p = ", format(p_value, scientific = TRUE, digits = 3)), "NS")
    # p_value_text <- ifelse(p_value_kw < 0.05, paste0("p = ", format(p_value_kw, scientific = TRUE, digits = 3)), "NS")
  }
  
  
  
  p1 <- ggplot(aes(x = filename, y = mean), data = current_cond) +
    geom_violin(draw_quantiles = c(0.25, 0.50, 0.75)) +
    geom_sina(aes(colour = filename))+
    theme(
      panel.background = element_blank(),
      axis.title = element_text(size = 18),
      axis.text = element_text(size = 16),
      axis.ticks.length = unit(0.2, units = "cm"),
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
      axis.line = element_line(color = "black", lineend = "round", linewidth = 1.25),
      legend.title = element_blank(),
      legend.text = element_blank(),
      legend.key = element_blank(),
      legend.position = "none"
    ) +
    ggtitle(paste0("K = 5 | ", current_cond$condition, " | ", p_value_text)) +
    xlab("Image") +
    ylab("Mean Intensity of Cells") +
    coord_flip()
  all_imgs[[counter]] <- p1
  counter <- counter + 1
}



p2 <- ggarrange(all_imgs[[4]],
                all_imgs[[3]],
                all_imgs[[2]],
                all_imgs[[1]],
                ncol = 2,
                nrow = 2,
                align = "hv",
                labels = "AUTO"
)


# Now doing scatter plots
plots <- list()
counter <- 1

for (c in unique(mean_neighbor_df_simplified$condition)) {
  current_cond <- filter(mean_neighbor_df_simplified, condition == c)
  current_cond$filename <- as.factor(current_cond$filename)
  p1 <- ggscatter(
    data = current_cond,
    x = "mean", y = "positive_neighborhood_mean",
    conf.int = TRUE,
    cor.coef = TRUE,
    cor.coeff.args = list(
      method = "spearman",
      label.x = 1000, label.y = 3000, label.sep = ",",
      cor.coef.name = "rho"
    ),
    parse = TRUE,
    cor.coef.size = 10,
    shape = 21,
    size = 2.5,
    add = fit_line,
    add.params = list(
      color = "black",
      fill = "grey"
    ),
    repel = TRUE,
    label.rectangle = TRUE,
    color = "filename",
    fill = "filename",
    alpha = 1,
    cor.method = stat_used,
    xlab = "Mean Intensity of Cells",
    ylab = "Mean MECP2+\nNeighborhood Intensity",
    ellipse = FALSE,
    mean.point = FALSE
  )
  
  # p2 for the single annotation context
  p2 <- ggpar(p1,
              title = paste0("K = ", k_n, " | ", c, " | 3D "),
              font.main = c(20, "bold"),
              font.x = c(18, "bold"),
              font.y = c(18, "bold"),
              ggtheme = theme(plot.title = element_text(hjust = 0.5))
  )
  
  p2 <- p2 +
    theme(
      axis.line = element_line(linewidth = 1.25, lineend = "round"),
      axis.text = element_text(size = 18),
      axis.ticks = element_line(linewidth = 0.2),
      plot.subtitle = element_text(face = "italic"),
      legend.title = element_text(size = 16, face = "bold"),
      legend.text = element_text(size = 14),
      legend.position = "bottom"
    )+
    labs(color = "Image", fill = "Image")
  
  # Change to p2 for single annotation case, p5 otherwise
  plots[[paste0(c)]] <- p2
  counter <- counter + 1
}

# Combo plot
combo_plot <- ggarrange(
  plots$NW,
  plots$NH,
  plots$SW,
  plots$SH,
  ncol = 2,
  nrow = 2,
  labels = "AUTO",
  align = "hv",
  legend = "right",
  common.legend = FALSE,
  font.label = list(
    size = 18,
    face = "bold"
  )
)
