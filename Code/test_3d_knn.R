# ----Loading packages----
pacman::p_load(
  car, # Tools for performing specific statistical tests
  caret, # Tools for training and evaluating predictive models
  class, # Functions for nearest neighbor classification
  deldir, # Delaunay triangulation and Dirichlet (Voronoi) tesselation
  doParallel, # Additional parallel functions
  effectsize, # Functions for calculating standardized effect sizes
  foreach, # Functions for parallel computing
  ggforce, # Additional geometries, stats, scales and themes for ggplot2
  ggpubr, # Functions for combining multiple ggplots into a single plot
  ggsignif, # Geometries for adding significance bars to ggplot2 plots
  parallel, # More functions for parallel computing
  progress, # Functions for progress bar
  rstatix, # Functions for descriptive statistics and statistical tests
  readxl, # Functions for reading Excel files
  reticulate, # Interface to Python for calling Python code from R
  scales, # Functions for scaling and formatting plot axes
  tidyverse # Collection of packages for data manipulation and visualization
)

# ----Loading 6 wk data----
six_wk_data <- read.csv("Data/Old data/Jacob/3D_MECP2_6WO_combined_dataset.csv")
six_wk_data <- six_wk_data %>%
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

six_wk_data <- six_wk_data %>% group_by(filename)
all_imgs <- six_wk_data %>% group_split(six_wk_data)

# ----Loading plotting function made in improved_knn_permutation.R----
#' Plot histogram for permutation test results
#'
#' This function creates a histogram plot based on the provided data, showing the results of a permutation test.
#'
#' @param data A data frame containing the permutation test results.
#' @param filename The base filename for saving the plot (without extension).
#' @param x_intercept The x-coordinate of the vertical red line to indicate the observed value.
#' @param plot_size The size of the plot (in inches). Default is 10.
#' @param base_plot_size The base size of the plot for calculating text scaling. Default is 8.
#' @param dpi The resolution for saving the plot (dots per inch). Default is 600.
#' @param hist_color The color of the histogram border. Default is "steelblue".
#' @param hist_fill The fill color of the histogram. Default is "steelblue".
#' @param num_bins The number of bins in the histogram. Default is 30.
#' @param line_color The color of the vertical red line. Default is "red".
#' @param text_color The color of the p-value text. Default is "red".
#' @param calculate_pvalue Logical indicating whether to calculate the p-value. Default is TRUE.
#' @param plot_pvalue Logical indicating whether to plot the p-value on the histogram. Default is FALSE.
#'
#' @return The created ggplot object representing the histogram plot.
#'
#' @import ggplot2
#' @import grid
#' @import scales
#' @import ggrepel
#'
#' @export
plot_permutation_test <- function(data,
                                  filename,
                                  x_intercept,
                                  plot_size = 10,
                                  base_plot_size = 8,
                                  dpi = 600,
                                  hist_color = "steelblue",
                                  hist_fill = "steelblue",
                                  num_bins = 30,
                                  line_color = "red",
                                  text_color = "red",
                                  calculate_pvalue = TRUE,
                                  p_value_test = "less",
                                  test_type = "One Sample t-test",
                                  test_direction = "Lesser") {
  # Required packages
  library(ggplot2)
  library(grid)
  library(scales)

  # Extracting the relevant pieces of data that will be used in building the plot
  time_point <- unique(data$time_point)
  k_value <- unique(data$k_value)
  num_permutations <- unique(data$num_iterations)
  num_ks <- unique(data$search_space)

  # Modifying the title of plots that are generated with max randomness
  if (num_ks == "max_randomness") {
    num_ks <- "All"
  }


  # Calculate the scaling factor for text elements based on plot size and base plot size
  text_scale <- plot_size / base_plot_size

  # Calculate p-value if required
  if (calculate_pvalue) {
    # Get the observed value (the x_intercept)
    observed_value <- x_intercept

    # Calculate the p-value
    p_value_comp <- t.test(data$rho, mu = observed_value, alternative = p_value_test)

    # Print the observed value and p-value to the console
    message(paste0("Observed Value: ", round(observed_value, digits = 4)))
    print(p_value_comp)

    # Build the plot with p-value if calculation is requested
    p <- ggplot(data = data, aes(x = rho)) +
      geom_histogram(color = hist_color, fill = hist_fill, bins = num_bins) +
      theme(
        panel.background = element_blank(),
        axis.title = element_text(size = 18 * text_scale, face = "bold"),
        axis.text = element_text(size = 16 * text_scale),
        axis.ticks = element_line(linewidth = 1.25 * text_scale),
        axis.ticks.length = unit(0.2, units = "cm") * text_scale,
        strip.text = element_text(size = 18 * text_scale, face = "bold"),
        plot.title = element_text(size = 20 * text_scale, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 18 * text_scale, face = "bold", hjust = 0.5),
        axis.line = element_line(colour = "black", linewidth = 1.25 * text_scale)
      ) +
      ggtitle(paste0("K = ", k_value, " | ", unique(data$condition), " | ", num_permutations, " Permutations | \n", time_point, " WK | ", num_ks, " KNNs Searched")) +
      xlab(expression(rho)) +
      ylab("Count") +
      scale_y_continuous(expand = c(0, 0)) +
      geom_vline(xintercept = x_intercept, color = line_color) +
      labs(subtitle = paste(test_type, " (", test_direction, ")\np-value: ", round(p_value_comp$p.value, digits = 4), sep = ""))
  } else {
    # Without p-value if not requested
    p <- ggplot(data = data, aes(x = rho)) +
      geom_histogram(color = hist_color, fill = hist_fill, bins = num_bins) +
      theme(
        panel.background = element_blank(),
        axis.title = element_text(size = 18 * text_scale, face = "bold"),
        axis.text = element_text(size = 16 * text_scale),
        axis.ticks = element_line(linewidth = 1.25 * text_scale),
        axis.ticks.length = unit(0.2, units = "cm") * text_scale,
        strip.text = element_text(size = 18 * text_scale, face = "bold"),
        plot.title = element_text(size = 20 * text_scale, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 18 * text_scale, face = "bold", hjust = 0.5),
        axis.line = element_line(colour = "black", linewidth = 1.25 * text_scale)
      ) +
      ggtitle(paste0("K = ", k_value, " | ", unique(data$condition), " | ", num_permutations, " Permutations | \n", time_point, " WK | ", num_ks, " KNNs Searched")) +
      xlab(expression(rho)) +
      ylab("Count") +
      scale_y_continuous(expand = c(0, 0)) +
      geom_vline(xintercept = x_intercept, color = line_color) +
      labs(subtitle = paste("No T-Test Performed", sep = ""))
  }


  # Calculate the width and height for the plot
  width <- plot_size * dpi
  height <- plot_size * dpi

  # If calculated p-value is less than 2.2e-16 and we are calculating the p-value than add subtitle mentioning raw p-value
  if (calculate_pvalue) {
    if (p_value_comp$p.value < 2.2e-16) {
      p <- p + labs(subtitle = paste(test_type, " (", test_direction, ")\np-value: ", round(p_value_comp$p.value, digits = 4), "\nRaw p-value < 2.2e-16", sep = ""))
    }
  }

  # Convert width and height to inches
  width_inches <- width / dpi
  height_inches <- height / dpi

  # Print width and height to the console
  message(paste0("Saving plot ", filename, " with width equal to ", width_inches, " inches and height equal to ", height_inches, " inches"))

  # Save as SVG
  ggsave(plot = p, filename = paste0(filename, ".svg"), device = "svg", width = width, height = height, units = "px", dpi = dpi, path = "Outputs/knn_analysis/plots/")

  # Save as PNG
  ggsave(plot = p, filename = paste0(filename, ".png"), device = "png", width = width, height = height, units = "px", dpi = dpi, path = "Outputs/knn_analysis/plots/")

  # Return the plot
  return(p)
}
# ----Function to calculate Euclidean distance between two 3D points----
euclidean_distance_3d <- function(x1, y1, z1, x2, y2, z2) {
  sqrt((x2 - x1)^2 + (y2 - y1)^2 + (z2 - z1)^2)
}

# ----Semi-random parallel code----
# Function to process one iteration
process_iteration <- function(iteration, max_iterations = 1000, time_point_to_use = 6, all_imgs, num_nearest, k_neighbors = 5, apply_intensity = TRUE, p_only = TRUE, calc_mean_intensity = TRUE, conditions = c("NW", "SW", "NH", "SH")) {
  
  random_seed <- runif(1, max_iterations, 100000)
  set.seed(random_seed)
  
  correlation_values <- list()
  counter <- 1
  iteration_data <- list()

  
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

      if (length(distances) >= num_nearest) {
        num_nearest_indices <- order(distances)[1:num_nearest]
        num_nearest_points <- other_points[num_nearest_indices, ]
      } else {
        next
      }

      if (num_nearest == k_neighbors) {
        sampled_indices <- sample(1:nrow(num_nearest_points), k_neighbors, replace = TRUE)
        k_nearest_points <- num_nearest_points[sampled_indices,]
      } else {
        sampled_indices <- sample(1:nrow(num_nearest_points), k_neighbors, replace = TRUE)
        k_nearest_points <- num_nearest_points[sampled_indices, ]
      }

      intensity_info <- k_nearest_points
      add_p_score <- sum((intensity_info$mecp2_p == "P")) / k_neighbors
      observation$knn_5_label <- add_p_score
      neighbor_means <- c()

      if (apply_intensity) {
        distances <- distances[num_nearest_indices]
        distances_updated <- scales::rescale(distances * intensity_info$mean, to = c(0, 1))
        distances_updated <- sort(distances_updated, decreasing = FALSE)

        if (p_only) {
          intensity_info <- dplyr::filter(intensity_info, mecp2_p == "P")
          if (nrow(intensity_info) == 0) {
            writeLines("This sample has only N neighbors.\nIt is being excluded from this analysis.")
          }
        }

        if (calc_mean_intensity) {
          df$positive_neighborhood_mean <- rep(0, nrow(df))
          mean_of_neighbors <- mean(intensity_info$mean)
          neighbor_means <- c(neighbor_means, mean_of_neighbors)
          observation$positive_neighborhood_mean <- mean_of_neighbors
          iteration_data[[counter]] <- observation
          counter <- counter + 1
        }
      }
    }
  }

  combined_data <- dplyr::bind_rows(iteration_data)

  # Set up default values
  for (cond in conditions) {
    correlation_values[[cond]] <- NA
  }

  conditions <- unique(combined_data$condition)

  for (cond in conditions) {
    specific_data <- dplyr::filter(combined_data, .data$condition == cond)
    valid_observations_mean <- sum(!is.na(specific_data$mean) & is.finite(specific_data$mean))
    valid_observations_positive <- sum(!is.na(specific_data$positive_neighborhood_mean) & is.finite(specific_data$positive_neighborhood_mean))


    if (valid_observations_mean > 1 && valid_observations_positive > 1) {
      correlation_test <- cor.test(specific_data$mean, specific_data$positive_neighborhood_mean, method = "spearman", use = "complete.obs")
      combined_data$search_space <- rep(num_nearest, nrow(combined_data))
      combined_data$k_value <- rep(k_neighbors, nrow(combined_data))
      combined_data$time_point <- rep(6, nrow(combined_data))
      combined_data$num_iterations <- rep(max_iterations, nrow(combined_data))
      
      current_df <- data.frame(
        search_space = rep(num_nearest, nrow(specific_data)),
        k_value = rep(k_neighbors, nrow(specific_data)),
        time_point = rep(time_point_to_use, nrow(specific_data)),
        num_iterations = rep(max_iterations, nrow(specific_data)),
        rho = rep(correlation_test$estimate, nrow(specific_data))
      )
      # correlation_values[[cond]] <- correlation_test$estimate
      correlation_values[[cond]] <- current_df
    }
  }

  # return(list(paste("Iteration", iteration) = correlation_values))
  return(correlation_values)
}


# Set the number of CPU cores to use (you can adjust this based on your system)
num_cores <- 11

# Initialize parallel backend
cl <- makeCluster(num_cores)
registerDoParallel(cl)

# Parallelized loop using foreach for an individual search space value
# total_iterations <- 4
# all_results <- list()
# counter <- 1
# search_space <- seq(5,50,5)
# for (s in search_space) {
#   for (i in 1:1000) {
#     results <- foreach(iteration = 1:total_iterations, .combine = c) %dopar% {
#       process_iteration(iteration, max_iterations = total_iterations, all_imgs = all_imgs, num_nearest = s, k_neighbors = 5, apply_intensity = TRUE, p_only = TRUE, calc_mean_intensity = TRUE, conditions = TRUE)
#     }
#     
#     all_results[[counter]] <- results
#     counter <- counter + 1
#   }
# }
# 
# 
# # Close the parallel backend
# stopCluster(cl)


# With progress bar
total_iterations <- 1
all_results <- list()
counter <- 1
search_space <- seq(5,50,5)

# Define the total number of tasks
total_tasks <- length(search_space) * 1000

# Create a progress bar
pb <- txtProgressBar(min = 0, max = total_tasks, style = 3)

for (s in search_space) {
  for (i in 1:1000) {
    results <- foreach(iteration = 1:total_iterations, .combine = c) %dopar% {
      process_iteration(iteration, max_iterations = total_iterations, all_imgs = all_imgs, num_nearest = s, k_neighbors = 5, apply_intensity = TRUE, p_only = TRUE, calc_mean_intensity = TRUE, conditions = TRUE)
    }
    
    # Remove elements with empty names
    results <- results[names(results) != ""]
    
    all_results[[paste0("Iteration ", counter, " Search Space ", s)]] <- results
    counter <- counter + 1
    
    # Update the progress bar
    setTxtProgressBar(pb, counter)
  }
}

# Close the progress bar
close(pb)

# Close the parallel backend
stopCluster(cl)




# x-intercepts (6 wk): NW: 0.27, NH: 0.44, SW: 0.52, SH: 0.58
conds <- c("NW", "NH", "SW", "SH")
x_ints <- c(0.27, 0.44, 0.52, 0.58)

# ----Plotting semi-random 1000 rep results----
for (cond in conds) {
  for (x in x_ints) {
    plot_permutation_test(data = results[[paste0(cond)]], filename = paste0(c,"_test2_50_nn"), x_intercept = x, calculate_pvalue = FALSE)
  }
}



# ----Completely random----
total_iterations <- 1000
k_neighbors <- 5
num_nearest <- 20
correlation_values <- list()

# Setting up the progress bar
pb <- txtProgressBar(min = 0, max = total_iterations, style = 3)

# Begin total_iterations
for (iteration in 1:total_iterations) {
  iteration_data <- list() # Data for this iteration
  counter <- 1

  for (i in 1:length(all_imgs)) {
    # Going through each image
    df <- all_imgs[[i]]
    df$knn_5_label <- rep(NA, nrow(df)) # Initialize the column
    df$positive_neighborhood_mean <- rep(NA, nrow(df)) # Initialize the column

    # For every row / observation within an image
    for (r in 1:nrow(df)) {
      observation <- df[r, ]
      other_points <- df[-r, ]

      distances <- euclidean_distance_3d(
        observation$x, observation$y, observation$z,
        other_points$x, other_points$y, other_points$z
      )

      # Randomly selecting from all the other points
      random_indices <- sample(1:nrow(other_points), k_neighbors)
      random_points <- other_points[random_indices, ]

      intensity_info <- random_points

      add_p_score <- sum((intensity_info$mecp2_p == "P")) / k_neighbors
      df$knn_5_label[r] <- add_p_score

      # Calculate distances with intensity information
      if (apply_intensity) {
        distances <- distances[random_indices]
        distances_updated <- rescale(distances * intensity_info$mean, to = c(0, 1))
        distances_updated <- sort(distances_updated, decreasing = FALSE)

        if (p_only) {
          intensity_info <- filter(intensity_info, mecp2_p == "P")
        }
        if (calc_mean_intensity) {
          mean_of_neighbors <- mean(intensity_info$mean)
          df$positive_neighborhood_mean[r] <- mean_of_neighbors
        }
      }
    }

    iteration_data[[i]] <- df
  }

  # Collate data and compute correlations for each condition
  combined_data <- bind_rows(iteration_data)
  conditions <- unique(combined_data$condition)

  correlation_results <- list()
  for (cond in conditions) {
    nw_data <- filter(combined_data, condition == cond)
    correlation_test <- cor.test(nw_data$mean, nw_data$positive_neighborhood_mean, method = "spearman", use = "complete.obs")
    correlation_results[[cond]] <- correlation_test$estimate
  }

  correlation_values[[paste("Iteration", iteration)]] <- correlation_results

  # Update the progress bar
  setTxtProgressBar(pb, iteration)
}

# Close the progress bar
close(pb)


# Filtering to just conditions to take a look at histograms
# Extracting all rhos for each iteration and condition
all_rhos <- lapply(correlation_values, function(iteration_result) {
  lapply(iteration_result, function(cond_result) {
    cond_result["rho"]
  })
})


# Convert list of lists to a matrix
rho_matrix <- do.call(rbind, lapply(all_rhos, unlist))

# SH
hist(rho_matrix[, 1], col = "lightblue", main = "6 WK SH", xlab = "Rho", xlim = c(0.38, 0.60))
abline(v = 0.58, col = "red", lwd = 3, lty = 1)


# SW
hist(rho_matrix[, 2], col = "lightblue", main = "6 WK SW", xlab = "Rho", xlim = c(0.33, 0.53))
abline(v = 0.52, col = "red", lwd = 3, lty = 1)


# NH
hist(rho_matrix[, 3], col = "lightblue", main = "6 WK NH", xlab = "Rho", xlim = c(0.18, 0.45))
abline(v = 0.44, col = "red", lwd = 3, lty = 1)

# NW
hist(rho_matrix[, 4], col = "lightblue", main = "6 WK NW", xlab = "Rho", xlim = c(-0.04, 0.28))
abline(v = 0.27, col = "red", lwd = 3, lty = 1)

# + Neighborhood mean
test_performed <- "kw"
plot_type <- "svg"
time_point <- "6wk"
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
    if (test_performed == "fk") {
      p_value_text <- ifelse(p_value < 0.05, paste0("p = ", format(p_value, scientific = TRUE, digits = 3)), "NS")
    } else {
      p_value_text <- ifelse(p_value_kw < 0.05, paste0("p = ", format(p_value_kw, scientific = TRUE, digits = 3)), "NS")
    }
  }



  p1 <- ggplot(aes(x = filename, y = positive_neighborhood_mean), data = current_cond) +
    geom_violin(draw_quantiles = c(0.25, 0.50, 0.75)) +
    geom_sina(aes(colour = filename)) +
    theme(
      panel.background = element_blank(),
      axis.title = element_text(size = 18, face = "bold"),
      axis.text = element_text(size = 16),
      axis.ticks.length = unit(0.2, units = "cm"),
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
      axis.line = element_line(color = "black", lineend = "round", linewidth = 1.25),
      legend.title = element_blank(),
      legend.text = element_blank(),
      legend.key = element_blank(),
      legend.position = "none"
    ) +
    ggtitle(paste0("K = 5 | ", current_cond$condition, " | ", toupper(test_performed), " | ", p_value_text, " | ", toupper(time_point))) +
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
  labels = "AUTO", font.label = list(size = 18, face = "bold")
)


# Saving combo plot
ggsave(p2,
  device = plot_type, path = "Outputs/",
  filename = paste0("positive_neighborhood_mean_violin_plot_combo_plot_3d_", test_performed, "_", time_point, ".", plot_type),
  width = 10, height = 10, units = "in", dpi = 600, bg = "white"
)


# Saving each individual plot
for (i in 1:length(all_imgs)) {
  ggsave(
    plot = all_imgs[[i]],
    filename = paste0("positive_neighborhood_mean_violin_plot_3d_", test_performed, "_", time_point, ".", plot_type),
    width = 10, height = 10, units = "in", dpi = 600, bg = "white", device = plot_type, path = "Outputs/"
  )
}


# Mean Intensity of Cells
test_performed <- "kw"
plot_type <- "svg"
time_point <- "6wk"
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
    if (test_performed == "fk") {
      p_value_text <- ifelse(p_value < 0.05, paste0("p = ", format(p_value, scientific = TRUE, digits = 3)), "NS")
    } else {
      p_value_text <- ifelse(p_value_kw < 0.05, paste0("p = ", format(p_value_kw, scientific = TRUE, digits = 3)), "NS")
    }
  }



  p1 <- ggplot(aes(x = filename, y = mean), data = current_cond) +
    geom_violin(draw_quantiles = c(0.25, 0.50, 0.75)) +
    geom_sina(aes(colour = filename)) +
    theme(
      panel.background = element_blank(),
      axis.title = element_text(size = 18, face = "bold"),
      axis.text = element_text(size = 16),
      axis.ticks.length = unit(0.2, units = "cm"),
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
      axis.line = element_line(color = "black", lineend = "round", linewidth = 1.25),
      legend.title = element_blank(),
      legend.text = element_blank(),
      legend.key = element_blank(),
      legend.position = "none"
    ) +
    ggtitle(paste0("K = 5 | ", current_cond$condition, " | ", toupper(test_performed), " | ", p_value_text, " | ", toupper(time_point))) +
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
  labels = "AUTO", font.label = list(size = 18, face = "bold")
)


# Saving combo plot
ggsave(p2,
  device = plot_type, path = "Outputs/",
  filename = paste0("mean_intensity_of_cells_violin_plot_combo_plot_3d_", test_performed, "_", time_point, ".", plot_type),
  width = 10, height = 10, units = "in", dpi = 600, bg = "white"
)


# Saving each individual plot
for (i in 1:length(all_imgs)) {
  ggsave(
    plot = all_imgs[[i]],
    filename = paste0("mean_intensity_of_cells_violin_plot_3d_", test_performed, "_", time_point, ".", plot_type),
    width = 10, height = 10, units = "in", dpi = 600, bg = "white", device = plot_type, path = "Outputs/"
  )
}


# Now doing scatter plots
plot_type <- "png"
time_point <- "6wk"
fit_line <- "reg.line"
stat_used <- "spearman"
plots <- list()
k_n <- 5
counter <- 1

mean_neighbor_df_simplified <- mean_neighbor_df_simplified %>%
  rename(Image = filename)

for (c in unique(mean_neighbor_df_simplified$condition)) {
  current_cond <- filter(mean_neighbor_df_simplified, condition == c)
  current_cond$Image <- as.factor(current_cond$Image)
  current_cond <- filter(current_cond, mecp2_p == "P")
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
    color = "Image",
    fill = "Image",
    alpha = 1,
    cor.method = stat_used,
    xlab = "Mean Intensity of Cells",
    ylab = "Mean MECP2+\nNeighborhood Intensity",
    ellipse = FALSE,
    mean.point = FALSE
  )

  # p2 for the single annotation context
  p2 <- ggpar(p1,
    title = paste0("K = ", k_n, " | ", c, " | 3D | ", toupper(time_point)),
    font.main = c(18, "bold"),
    font.x = c(18, "bold"),
    font.y = c(18, "bold"),
    ggtheme = theme(plot.title = element_text(hjust = 0.5))
  )

  p2 <- p2 +
    theme(
      axis.line = element_line(linewidth = 1.25, lineend = "round"),
      axis.title = element_text(size = 18, face = "bold"),
      axis.text = element_text(size = 16),
      axis.ticks = element_line(linewidth = 0.2),
      plot.subtitle = element_text(face = "italic"),
      legend.title = element_text(size = 16, face = "bold"),
      legend.text = element_text(size = 14),
      legend.position = "bottom"
    ) +
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


# Saving combo plot
ggsave(combo_plot,
  device = plot_type, path = "Outputs/",
  filename = paste0("scatter_plot_3d_combo_", "_", time_point, ".", plot_type),
  width = 14, height = 14, units = "in", dpi = 600, bg = "white"
)


# Saving each individual plot
for (i in 1:length(plots)) {
  ggsave(
    plot = plots[[i]],
    filename = paste0("scatter_plot_3d_", i, "_", time_point, ".", plot_type),
    width = 10, height = 10, units = "in", dpi = 600, bg = "white", device = plot_type, path = "Outputs/"
  )
}



# 12 week dataset
twelve_df <- read.csv("Data/Old data/Jacob/12WOcombined_output_MECP2.csv")
twelve_df <- twelve_df %>%
  select(Image, CX..pix., CY..pix., CZ..pix., Condition, Hemishphere, MECP2, Mean) %>%
  rename(
    x = CX..pix.,
    y = CY..pix.,
    z = CZ..pix.,
    mecp2_p = MECP2
  ) %>%
  rename_all(tolower) %>%
  filter(mecp2_p == "P")

twelve_df <- twelve_df %>% group_by(image)
all_imgs <- twelve_df %>% group_split(twelve_df)


# KNN analysis
all_data <- list()
all_p_scores <- list()
counter <- 1
apply_intensity <- TRUE
p_only <- TRUE
calc_mean_intensity <- TRUE
neighbor_means <- c()


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
# test_df <- test_df %>%
#   rename(filename = image)

# + Neighborhood mean
test_performed <- "kw"
plot_type <- "svg"
time_point <- "12wk"
all_imgs <- list()
counter <- 1
mean_neighbor_df_simplified <- test_df
mean_neighbor_df_simplified <- mean_neighbor_df_simplified %>%
  mutate(image = match(image, unique(image)))

for (i in unique(mean_neighbor_df_simplified$condition)) {
  current_cond <- filter(mean_neighbor_df_simplified, condition == i)
  current_cond <- filter(current_cond, mecp2_p == "P")
  current_cond$image <- as.factor(current_cond$image)

  # Test for normality using Shapiro-Wilk test
  shapiro_test <- shapiro.test(current_cond$positive_neighborhood_mean)
  if (shapiro_test$p.value >= 0.05) {
    # Data is normally distributed
    print("Data is normally distributed")

    # Perform Levene's test for homogeneity of variances
    levene_test <- leveneTest(positive_neighborhood_mean ~ image, data = current_cond)
    print(levene_test)

    p_value <- levene_test$"Pr(>F)"[1] # Get the p-value from the test
    p_value_text <- ifelse(p_value < 0.05, paste0("p = ", format(p_value, scientific = FALSE, digits = 3)), "NS")
  } else {
    # Data is not normally distributed
    print("Data is not normally distributed")

    # Perform Fligner-Killeen test for homogeneity of variances
    fligner_test <- fligner.test(positive_neighborhood_mean ~ image, data = current_cond)
    print(fligner_test)


    kw_test <- kruskal.test(positive_neighborhood_mean ~ image, data = current_cond)
    print(kw_test)

    p_value <- fligner_test$p.value
    p_value_kw <- kw_test$p.value
    if (test_performed == "fk") {
      p_value_text <- ifelse(p_value < 0.05, paste0("p = ", format(p_value, scientific = TRUE, digits = 3)), "NS")
    } else {
      p_value_text <- ifelse(p_value_kw < 0.05, paste0("p = ", format(p_value_kw, scientific = TRUE, digits = 3)), "NS")
    }
  }



  p1 <- ggplot(aes(x = image, y = positive_neighborhood_mean), data = current_cond) +
    geom_violin(draw_quantiles = c(0.25, 0.50, 0.75)) +
    geom_sina(aes(colour = image)) +
    theme(
      panel.background = element_blank(),
      axis.title = element_text(size = 18, face = "bold"),
      axis.text = element_text(size = 16),
      axis.ticks.length = unit(0.2, units = "cm"),
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
      axis.line = element_line(color = "black", lineend = "round", linewidth = 1.25),
      legend.title = element_blank(),
      legend.text = element_blank(),
      legend.key = element_blank(),
      legend.position = "none"
    ) +
    ggtitle(paste0("K = 5 | ", current_cond$condition, " | ", toupper(test_performed), " | ", p_value_text, " | ", toupper(time_point))) +
    xlab("Image") +
    ylab("+ Neighborhood Mean") +
    coord_flip()
  all_imgs[[counter]] <- p1
  counter <- counter + 1
}



p2 <- ggarrange(all_imgs[[2]],
  all_imgs[[1]],
  all_imgs[[4]],
  all_imgs[[3]],
  ncol = 2,
  nrow = 2,
  align = "hv",
  labels = "AUTO", font.label = list(size = 18, face = "bold")
)


# Saving combo plot
ggsave(p2,
  device = plot_type, path = "Outputs/",
  filename = paste0("positive_neighborhood_mean_violin_plot_combo_plot_3d_", test_performed, "_", time_point, ".", plot_type),
  width = 10, height = 10, units = "in", dpi = 600, bg = "white"
)


# Saving each individual plot
for (i in 1:length(all_imgs)) {
  ggsave(
    plot = all_imgs[[i]],
    filename = paste0("positive_neighborhood_mean_violin_plot_3d_", test_performed, "_", time_point, ".", plot_type),
    width = 10, height = 10, units = "in", dpi = 600, bg = "white", device = plot_type, path = "Outputs/"
  )
}


# Mean Intensity of Cells
test_performed <- "kw"
plot_type <- "png"
time_point <- "12wk"
all_imgs <- list()
counter <- 1
mean_neighbor_df_simplified <- test_df
mean_neighbor_df_simplified <- mean_neighbor_df_simplified %>%
  mutate(image = match(image, unique(image)))
mean_neighbor_df_simplified <- filter(mean_neighbor_df_simplified, mecp2_p == "P")

for (i in unique(mean_neighbor_df_simplified$condition)) {
  current_cond <- filter(mean_neighbor_df_simplified, condition == i)
  current_cond$image <- as.factor(current_cond$image)

  # Test for normality using Shapiro-Wilk test
  shapiro_test <- shapiro.test(current_cond$mean)
  if (shapiro_test$p.value >= 0.05) {
    # Data is normally distributed
    print("Data is normally distributed")

    # Perform Levene's test for homogeneity of variances
    levene_test <- leveneTest(mean ~ image, data = current_cond)
    print(levene_test)

    p_value <- levene_test$"Pr(>F)"[1] # Get the p-value from the test
    p_value_text <- ifelse(p_value < 0.05, paste0("p = ", format(p_value, scientific = FALSE, digits = 3)), "NS")
  } else {
    # Data is not normally distributed
    print("Data is not normally distributed")

    # Perform Fligner-Killeen test for homogeneity of variances
    fligner_test <- fligner.test(mean ~ image, data = current_cond)
    print(fligner_test)


    kw_test <- kruskal.test(mean ~ image, data = current_cond)
    print(kw_test)

    p_value <- fligner_test$p.value
    p_value_kw <- kw_test$p.value
    if (test_performed == "fk") {
      p_value_text <- ifelse(p_value < 0.05, paste0("p = ", format(p_value, scientific = TRUE, digits = 3)), "NS")
    } else {
      p_value_text <- ifelse(p_value_kw < 0.05, paste0("p = ", format(p_value_kw, scientific = TRUE, digits = 3)), "NS")
    }
  }



  p1 <- ggplot(aes(x = image, y = mean), data = current_cond) +
    geom_violin(draw_quantiles = c(0.25, 0.50, 0.75)) +
    geom_sina(aes(colour = image)) +
    theme(
      panel.background = element_blank(),
      axis.title = element_text(size = 18, face = "bold"),
      axis.text = element_text(size = 16),
      axis.ticks.length = unit(0.2, units = "cm"),
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
      axis.line = element_line(color = "black", lineend = "round", linewidth = 1.25),
      legend.title = element_blank(),
      legend.text = element_blank(),
      legend.key = element_blank(),
      legend.position = "none"
    ) +
    ggtitle(paste0("K = 5 | ", current_cond$condition, " | ", toupper(test_performed), " | ", p_value_text, " | ", toupper(time_point))) +
    xlab("Image") +
    ylab("Mean Intensity of Cells") +
    coord_flip()
  all_imgs[[counter]] <- p1
  counter <- counter + 1
}



p2 <- ggarrange(all_imgs[[2]],
  all_imgs[[1]],
  all_imgs[[4]],
  all_imgs[[3]],
  ncol = 2,
  nrow = 2,
  align = "hv",
  labels = "AUTO", font.label = list(size = 18, face = "bold")
)


# Saving combo plot
ggsave(p2,
  device = plot_type, path = "Outputs/",
  filename = paste0("mean_intensity_of_cells_violin_plot_combo_plot_3d_", test_performed, "_", time_point, ".", plot_type),
  width = 10, height = 10, units = "in", dpi = 600, bg = "white"
)


# Saving each individual plot
for (i in 1:length(all_imgs)) {
  ggsave(
    plot = all_imgs[[i]],
    filename = paste0("mean_intensity_of_cells_violin_plot_3d_", test_performed, "_", time_point, ".", plot_type),
    width = 10, height = 10, units = "in", dpi = 600, bg = "white", device = plot_type, path = "Outputs/"
  )
}


# Now doing scatter plots
plot_type <- "svg"
time_point <- "12wk"
plots <- list()
counter <- 1

mean_neighbor_df_simplified <- mean_neighbor_df_simplified %>%
  rename(Image = image)

for (c in unique(mean_neighbor_df_simplified$condition)) {
  current_cond <- filter(mean_neighbor_df_simplified, condition == c)
  current_cond$Image <- as.factor(current_cond$Image)
  current_cond <- filter(current_cond, mecp2_p == "P")
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
    color = "Image",
    fill = "Image",
    alpha = 1,
    cor.method = stat_used,
    xlab = "Mean Intensity of Cells",
    ylab = "Mean MECP2+\nNeighborhood Intensity",
    ellipse = FALSE,
    mean.point = FALSE
  )

  # p2 for the single annotation context
  p2 <- ggpar(p1,
    title = paste0("K = ", k_n, " | ", c, " | 3D | ", toupper(time_point)),
    font.main = c(18, "bold"),
    font.x = c(18, "bold"),
    font.y = c(18, "bold"),
    ggtheme = theme(plot.title = element_text(hjust = 0.5))
  )

  p2 <- p2 +
    theme(
      axis.line = element_line(linewidth = 1.25, lineend = "round"),
      axis.title = element_text(size = 18, face = "bold"),
      axis.text = element_text(size = 16),
      axis.ticks = element_line(linewidth = 0.2),
      plot.subtitle = element_text(face = "italic"),
      legend.title = element_text(size = 16, face = "bold"),
      legend.text = element_text(size = 14),
      legend.position = "bottom"
    ) +
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


# Saving combo plot
ggsave(combo_plot,
  device = plot_type, path = "Outputs/",
  filename = paste0("scatter_plot_3d_combo_", "_", time_point, ".", plot_type),
  width = 14, height = 14, units = "in", dpi = 600, bg = "white"
)


# Saving each individual plot
for (i in 1:length(plots)) {
  ggsave(
    plot = plots[[i]],
    filename = paste0("scatter_plot_3d_", i, "_", time_point, ".", plot_type),
    width = 10, height = 10, units = "in", dpi = 600, bg = "white", device = plot_type, path = "Outputs/"
  )
}


# Now looking to see if varying K from 1-15 changes for either 6 or 12 week
# 6 week data
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

# 12 week data
twelve_df <- read.csv("Data/Old data/Jacob/12WOcombined_output_MECP2.csv")
twelve_df <- twelve_df %>%
  select(Image, CX..pix., CY..pix., CZ..pix., Condition, Hemishphere, MECP2, Mean) %>%
  rename(
    x = CX..pix.,
    y = CY..pix.,
    z = CZ..pix.,
    mecp2_p = MECP2
  ) %>%
  rename_all(tolower) %>%
  filter(mecp2_p == "P")

twelve_df <- twelve_df %>% group_by(image)
all_imgs <- twelve_df %>% group_split(twelve_df)



ks <- seq(1, 15, 1)
time_point <- "12wk"
for (k in ks) {
  all_data <- list()
  all_p_scores <- list()
  counter <- 1
  apply_intensity <- TRUE
  p_only <- TRUE
  calc_mean_intensity <- TRUE
  neighbor_means <- c()


  # k_neighbors <- 5
  for (i in 1:length(all_imgs)) {
    df <- all_imgs[[i]]
    df$knn_label <- rep(0, nrow(df))
    for (r in 1:nrow(df)) {
      observation <- df[r, ]
      other_points <- df[-r, ]
      distances <- euclidean_distance_3d(
        observation$x, observation$y, observation$z,
        other_points$x, other_points$y, other_points$z
      )


      k_nearest_indices <- order(distances)[1:k]
      k_nearest_points <- other_points[k_nearest_indices, ]

      intensity_info <- other_points[k_nearest_indices, ]
      add_p_score <- sum((intensity_info$mecp2_p == "P")) / k
      observation$knn_label <- add_p_score
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
  write.csv(test_df, file = paste0("Outputs/", k, "_scatter_", time_point, ".csv"))
}

all_6wk_dfs <- list.files(path = "Outputs/", pattern = "*_6wk.csv")
all_12wk_dfs <- list.files(path = "Outputs/", pattern = "*_12wk.csv")

# 6 week rhos
all_rhos <- list()
counter <- 1
for (f in all_6wk_dfs) {
  current_k <- read.csv(paste0("Outputs/", f))
  for (c in unique(current_k$condition)) {
    current_k$rho <- rep(0, nrow(current_k))
    current_cond <- filter(current_k, condition == c)
    current_cond <- filter(current_cond, mecp2_p == "P")
    current_rho <- cor.test(x = current_cond$mean, y = current_cond$positive_neighborhood_mean, method = "spearman")
    all_rhos[[paste0(c, "_", f)]] <- current_rho$estimate
    counter <- counter + 1
  }
}


conditions <- unlist(str_extract_all(names(all_rhos), pattern = "(SH|SW|NH|NW)"))
rhos <- unlist(as.numeric(unname(all_rhos)))
ks <- unlist(str_extract_all(names(all_rhos), pattern = "\\d+"))
six_wk_rho_df <- data.frame(conditions, rhos, ks)
six_wk_rho_df <- six_wk_rho_df %>%
  arrange(., ks)

# 12 week rhos
all_rhos <- list()
counter <- 1
for (f in all_12wk_dfs) {
  current_k <- read.csv(paste0("Outputs/", f))
  for (c in unique(current_k$condition)) {
    current_k$rho <- rep(0, nrow(current_k))
    current_cond <- filter(current_k, condition == c)
    current_cond <- filter(current_cond, mecp2_p == "P")
    current_rho <- cor.test(x = current_cond$mean, y = current_cond$positive_neighborhood_mean, method = "spearman")
    all_rhos[[paste0(c, "_", f)]] <- current_rho$estimate
    counter <- counter + 1
  }
}


conditions <- unlist(str_extract_all(names(all_rhos), pattern = "(SH|SW|NH|NW)"))
rhos <- unlist(as.numeric(unname(all_rhos)))
ks <- unlist(str_extract_all(names(all_rhos), pattern = "\\d+"))
twelve_wk_rho_df <- data.frame(conditions, rhos, ks)
twelve_wk_rho_df <- twelve_wk_rho_df %>%
  arrange(., ks)


# 6 week plot
six_wk_rho_df$ks <- as.factor(six_wk_rho_df$ks)


p1 <- ggplot(data = six_wk_rho_df, aes(x = ks, y = rhos)) +
  geom_col(width = 0.8, fill = "steelblue") +
  theme(
    panel.background = element_blank(),
    axis.title = element_text(size = 18, face = "bold"),
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
    axis.line = element_line(color = "black", linewidth = 1.5),
    axis.ticks.length = unit(0.2, units = "cm"),
    axis.text = element_text(size = 16),
    legend.position = "bottom",
    legend.title = element_text(size = 16, face = "bold"),
    legend.text = element_text(size = 16),
    strip.text = element_text(size = 16, face = "bold"),
  ) +
  geom_text(aes(label = round(rhos, 2)), vjust = -0.5, size = 8) +
  facet_wrap(~conditions, ncol = 2, nrow = 2) +
  xlab("K") +
  ylab(expression(rho)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0.0, 0.62))

p1


ggsave(plot = p1, filename = "rhos_across_k_1_to_15_all_conds.svg", path = "Outputs/", device = "svg", height = 9, width = 16, units = "in", dpi = 600)

ggsave(plot = p1, filename = "rhos_across_k_1_to_15_all_conds.png", path = "Outputs/", device = "png", height = 9, width = 16, units = "in", dpi = 600)


# Permutation test for 3D 6 wk data
all_data <- list()
all_p_scores <- list()
counter <- 1
apply_intensity <- TRUE
p_only <- TRUE
calc_mean_intensity <- TRUE
neighbor_means <- c()
random_num_neighbors <- TRUE
k_neighbors <- 5
all_rhos_stored <- c()
all_p_values_stored <- c()


# 6 week data
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
  mutate(id = 1:nrow(.)) %>%
  # Change this line to filter to different condition
  filter(condition == "NW")

threed_data <- threed_data %>% group_by(filename)
all_imgs <- threed_data %>% group_split(threed_data)


for (n in 1:1000) {
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

      # Selecting random neighbors
      if (random_num_neighbors) {
        k_nearest_indices <- sample(nrow(other_points), k_neighbors, replace = TRUE)
        k_nearest_points <- other_points[k_nearest_indices, ]
      } else {
        distances <- euclidean_distance_3d(
          observation$x, observation$y, observation$z,
          other_points$x, other_points$y, other_points$z
        )
        k_nearest_indices <- order(distances)[1:k_neighbors]
        k_nearest_points <- other_points[k_nearest_indices, ]
      }

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
          df$positive_neighborhood_mean <- rep(0, nrow(df))
          df$rho <- rep(0, nrow(df))
          mean_of_neighbors <- mean(intensity_info$mean)
          observation$positive_neighborhood_mean <- mean_of_neighbors
          all_data[[counter]] <- observation
          counter <- counter + 1
        }
      }
    }
  }
  # Binding all observations from this group of images
  test_df <- bind_rows(all_data)

  # Calculating rho
  cor_res <- cor.test(x = test_df$mean, y = test_df$positive_neighborhood_mean, method = "spearman")
  current_p <- cor_res$p.value
  current_rho <- as.numeric(as.vector(unlist(cor_res$estimate)))
  all_rhos_stored <- c(all_rhos_stored, current_rho)
  all_p_values_stored <- c(all_p_values_stored, current_p)
}

nh_data <- data.frame(condition = rep("NH", times = length(all_rhos_stored)), rho = all_rhos_stored, p_value = all_p_values_stored)
nw_data <- data.frame(condition = rep("NW", times = length(all_rhos_stored)), rho = all_rhos_stored, p_value = all_p_values_stored)
sw_data <- data.frame(condition = rep("SW", times = length(all_rhos_stored)), rho = all_rhos_stored, p_value = all_p_values_stored)
sh_data <- data.frame(condition = rep("SH", times = length(all_rhos_stored)), rho = all_rhos_stored, p_value = all_p_values_stored)
# write.csv(nh_data, "Outputs/nh_1000_permutation_6wk_df.csv")
# write.csv(nw_data, "Outputs/nw_1000_permutation_6wk_df.csv")
# write.csv(sw_data, "Outputs/sw_1000_permutation_6wk_df.csv")
# write.csv(sh_data, "Outputs/sh_1000_permutation_6wk_df.csv")

# NH histogram
p1 <- ggplot(data = nh_data, aes(x = rho)) +
  geom_histogram(color = "steelblue", fill = "steelblue") +
  theme(
    panel.background = element_blank(),
    axis.title = element_text(size = 18, face = "bold"),
    axis.text = element_text(size = 16),
    axis.ticks = element_line(linewidth = 1.25),
    axis.ticks.length = unit(0.2, units = "cm"),
    strip.text = element_text(size = 18, face = "bold"),
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
    axis.line = element_line(colour = "black", linewidth = 1.25)
  ) +
  ggtitle(paste0("K = 5 | ", unique(nh_data$condition), " | ", length(all_rhos_stored))) +
  xlab(expression(rho)) +
  ylab("Count") +
  scale_y_continuous(expand = c(0, 0)) +
  geom_vline(xintercept = 0.44, color = "red")

# SVG
ggsave(plot = p1, filename = "nh_random_knn_5_vs_rho_histogram_6wk.svg", device = "svg", width = 10, height = 10, units = "in", path = "Outputs/knn_analysis/plots/")

# PNG
ggsave(plot = p1, filename = "nh_random_knn_5_vs_rho_histogram_6wk.png", device = "png", width = 10, height = 10, units = "in", path = "Outputs/knn_analysis/plots/")



# NW histogram
p1 <- ggplot(data = nw_data, aes(x = rho)) +
  geom_histogram(color = "steelblue", fill = "steelblue") +
  theme(
    panel.background = element_blank(),
    axis.title = element_text(size = 18, face = "bold"),
    axis.text = element_text(size = 16),
    axis.ticks = element_line(linewidth = 1.25),
    axis.ticks.length = unit(0.2, units = "cm"),
    strip.text = element_text(size = 18, face = "bold"),
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
    axis.line = element_line(colour = "black", linewidth = 1.25)
  ) +
  ggtitle(paste0("K = 5 | ", unique(nw_data$condition), " | ", length(all_rhos_stored))) +
  xlab(expression(rho)) +
  ylab("Count") +
  scale_y_continuous(expand = c(0, 0)) +
  geom_vline(xintercept = 0.27, color = "red")

# SVG
ggsave(plot = p1, filename = "nw_random_knn_5_vs_rho_histogram_6wk.svg", device = "svg", width = 10, height = 10, units = "in", path = "Outputs/knn_analysis/plots/")

# PNG
ggsave(plot = p1, filename = "nw_random_knn_5_vs_rho_histogram_6wk.png", device = "png", width = 10, height = 10, units = "in", path = "Outputs/knn_analysis/plots/")


# SW histogram
p1 <- ggplot(data = sw_data, aes(x = rho)) +
  geom_histogram(color = "steelblue", fill = "steelblue") +
  theme(
    panel.background = element_blank(),
    axis.title = element_text(size = 18, face = "bold"),
    axis.text = element_text(size = 16),
    axis.ticks = element_line(linewidth = 1.25),
    axis.ticks.length = unit(0.2, units = "cm"),
    strip.text = element_text(size = 18, face = "bold"),
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
    axis.line = element_line(colour = "black", linewidth = 1.25)
  ) +
  ggtitle(paste0("K = 5 | ", unique(sw_data$condition), " | ", length(all_rhos_stored))) +
  xlab(expression(rho)) +
  ylab("Count") +
  scale_y_continuous(expand = c(0, 0)) +
  geom_vline(xintercept = 0.52, color = "red")

# SVG
ggsave(plot = p1, filename = "sw_random_knn_5_vs_rho_histogram_6wk.svg", device = "svg", width = 10, height = 10, units = "in", path = "Outputs/knn_analysis/plots/")

# PNG
ggsave(plot = p1, filename = "sw_random_knn_5_vs_rho_histogram_6wk.png", device = "png", width = 10, height = 10, units = "in", path = "Outputs/knn_analysis/plots/")



# SH histogram
p1 <- ggplot(data = sh_data, aes(x = rho)) +
  geom_histogram(color = "steelblue", fill = "steelblue") +
  theme(
    panel.background = element_blank(),
    axis.title = element_text(size = 18, face = "bold"),
    axis.text = element_text(size = 16),
    axis.ticks = element_line(linewidth = 1.25),
    axis.ticks.length = unit(0.2, units = "cm"),
    strip.text = element_text(size = 18, face = "bold"),
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
    axis.line = element_line(colour = "black", linewidth = 1.25)
  ) +
  ggtitle(paste0("K = 5 | ", unique(sh_data$condition), " | ", length(all_rhos_stored))) +
  xlab(expression(rho)) +
  ylab("Count") +
  scale_y_continuous(expand = c(0, 0)) +
  geom_vline(xintercept = 0.58, color = "red")

# SVG
ggsave(plot = p1, filename = "sh_random_knn_5_vs_rho_histogram_6wk.svg", device = "svg", width = 10, height = 10, units = "in", path = "Outputs/knn_analysis/plots/")

# PNG
ggsave(plot = p1, filename = "sh_random_knn_5_vs_rho_histogram_6wk.png", device = "png", width = 10, height = 10, units = "in", path = "Outputs/knn_analysis/plots/")




# 12 week 3D data
all_data <- list()
all_p_scores <- list()
counter <- 1
apply_intensity <- TRUE
p_only <- TRUE
calc_mean_intensity <- TRUE
neighbor_means <- c()
random_num_neighbors <- TRUE
k_neighbors <- 5
all_rhos_stored <- c()
all_p_values_stored <- c()


twelve_df <- read.csv("Data/Old data/Jacob/12WOcombined_output_MECP2.csv")
twelve_df <- twelve_df %>%
  select(Image, CX..pix., CY..pix., CZ..pix., Condition, Hemishphere, MECP2, Mean) %>%
  rename(
    x = CX..pix.,
    y = CY..pix.,
    z = CZ..pix.,
    mecp2_p = MECP2
  ) %>%
  rename_all(tolower) %>%
  filter(mecp2_p == "P") %>%
  # Change this line to filter to other conditions
  filter(condition == "NW")

twelve_df <- twelve_df %>% group_by(image)
all_imgs <- twelve_df %>% group_split(twelve_df)
# Removing this image because it only contains 4 samples and can't be used for KNN = 5 analysis
# all_imgs <- all_imgs[-11]



for (n in 1:1000) {
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

      # Selecting random neighbors
      if (random_num_neighbors) {
        k_nearest_indices <- sample(nrow(other_points), k_neighbors, replace = TRUE)
        k_nearest_points <- other_points[k_nearest_indices, ]
      } else {
        distances <- euclidean_distance_3d(
          observation$x, observation$y, observation$z,
          other_points$x, other_points$y, other_points$z
        )
        k_nearest_indices <- order(distances)[1:k_neighbors]
        k_nearest_points <- other_points[k_nearest_indices, ]
      }

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
          df$positive_neighborhood_mean <- rep(0, nrow(df))
          df$rho <- rep(0, nrow(df))
          mean_of_neighbors <- mean(intensity_info$mean)
          observation$positive_neighborhood_mean <- mean_of_neighbors
          all_data[[counter]] <- observation
          counter <- counter + 1
        }
      }
    }
  }
  # Binding all observations from this group of images
  test_df <- bind_rows(all_data)

  # Calculating rho
  cor_res <- cor.test(x = test_df$mean, y = test_df$positive_neighborhood_mean, method = "spearman")
  current_p <- cor_res$p.value
  current_rho <- as.numeric(as.vector(unlist(cor_res$estimate)))
  all_rhos_stored <- c(all_rhos_stored, current_rho)
  all_p_values_stored <- c(all_p_values_stored, current_p)
}

nh_data <- data.frame(condition = rep("NH", times = length(all_rhos_stored)), rho = all_rhos_stored, p_value = all_p_values_stored)
nw_data <- data.frame(condition = rep("NW", times = length(all_rhos_stored)), rho = all_rhos_stored, p_value = all_p_values_stored)
sw_data <- data.frame(condition = rep("SW", times = length(all_rhos_stored)), rho = all_rhos_stored, p_value = all_p_values_stored)
sh_data <- data.frame(condition = rep("SH", times = length(all_rhos_stored)), rho = all_rhos_stored, p_value = all_p_values_stored)
# write.csv(nh_data, "Outputs/nh_1000_permutation_12wk_df.csv")
# write.csv(nw_data, "Outputs/nw_1000_permutation_12wk_df.csv")
# write.csv(sw_data, "Outputs/sw_1000_permutation_12wk_df.csv")
# write.csv(sh_data, "Outputs/sh_1000_permutation_12wk_df.csv")

# NH histogram
p1 <- ggplot(data = nh_data, aes(x = rho)) +
  geom_histogram(color = "steelblue", fill = "steelblue") +
  theme(
    panel.background = element_blank(),
    axis.title = element_text(size = 18, face = "bold"),
    axis.text = element_text(size = 16),
    axis.ticks = element_line(linewidth = 1.25),
    axis.ticks.length = unit(0.2, units = "cm"),
    strip.text = element_text(size = 18, face = "bold"),
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
    axis.line = element_line(colour = "black", linewidth = 1.25)
  ) +
  ggtitle(paste0("K = 5 | ", unique(nh_data$condition), " | ", length(all_rhos_stored), " | 12 WK")) +
  xlab(expression(rho)) +
  ylab("Count") +
  scale_y_continuous(expand = c(0, 0)) +
  geom_vline(xintercept = 0.17, color = "red")

# SVG
ggsave(plot = p1, filename = "nh_random_knn_5_vs_rho_histogram_12wk.svg", device = "svg", width = 10, height = 10, units = "in", path = "Outputs/knn_analysis/plots/")

# PNG
ggsave(plot = p1, filename = "nh_random_knn_5_vs_rho_histogram_12wk.png", device = "png", width = 10, height = 10, units = "in", path = "Outputs/knn_analysis/plots/")



# NW histogram
p1 <- ggplot(data = nw_data, aes(x = rho)) +
  geom_histogram(color = "steelblue", fill = "steelblue") +
  theme(
    panel.background = element_blank(),
    axis.title = element_text(size = 18, face = "bold"),
    axis.text = element_text(size = 16),
    axis.ticks = element_line(linewidth = 1.25),
    axis.ticks.length = unit(0.2, units = "cm"),
    strip.text = element_text(size = 18, face = "bold"),
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
    axis.line = element_line(colour = "black", linewidth = 1.25)
  ) +
  ggtitle(paste0("K = 5 | ", unique(nw_data$condition), " | ", length(all_rhos_stored), " | 12 WK")) +
  xlab(expression(rho)) +
  ylab("Count") +
  scale_y_continuous(expand = c(0, 0)) +
  geom_vline(xintercept = 0.27, color = "red")

# SVG
ggsave(plot = p1, filename = "nw_random_knn_5_vs_rho_histogram_12wk.svg", device = "svg", width = 10, height = 10, units = "in", path = "Outputs/knn_analysis/plots/")

# PNG
ggsave(plot = p1, filename = "nw_random_knn_5_vs_rho_histogram_12wk.png", device = "png", width = 10, height = 10, units = "in", path = "Outputs/knn_analysis/plots/")


# SW histogram
p1 <- ggplot(data = sw_data, aes(x = rho)) +
  geom_histogram(color = "steelblue", fill = "steelblue") +
  theme(
    panel.background = element_blank(),
    axis.title = element_text(size = 18, face = "bold"),
    axis.text = element_text(size = 16),
    axis.ticks = element_line(linewidth = 1.25),
    axis.ticks.length = unit(0.2, units = "cm"),
    strip.text = element_text(size = 18, face = "bold"),
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
    axis.line = element_line(colour = "black", linewidth = 1.25)
  ) +
  ggtitle(paste0("K = 5 | ", unique(sw_data$condition), " | ", length(all_rhos_stored), " | 12 WK")) +
  xlab(expression(rho)) +
  ylab("Count") +
  scale_y_continuous(expand = c(0, 0)) +
  geom_vline(xintercept = 0.52, color = "red")

# SVG
ggsave(plot = p1, filename = "sw_random_knn_5_vs_rho_histogram_12wk.svg", device = "svg", width = 10, height = 10, units = "in", path = "Outputs/knn_analysis/plots/")

# PNG
ggsave(plot = p1, filename = "sw_random_knn_5_vs_rho_histogram_12wk.png", device = "png", width = 10, height = 10, units = "in", path = "Outputs/knn_analysis/plots/")



# SH histogram
p1 <- ggplot(data = sh_data, aes(x = rho)) +
  geom_histogram(color = "steelblue", fill = "steelblue") +
  theme(
    panel.background = element_blank(),
    axis.title = element_text(size = 18, face = "bold"),
    axis.text = element_text(size = 16),
    axis.ticks = element_line(linewidth = 1.25),
    axis.ticks.length = unit(0.2, units = "cm"),
    strip.text = element_text(size = 18, face = "bold"),
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
    axis.line = element_line(colour = "black", linewidth = 1.25)
  ) +
  ggtitle(paste0("K = 5 | ", unique(sh_data$condition), " | ", length(all_rhos_stored), " | 12 WK")) +
  xlab(expression(rho)) +
  ylab("Count") +
  scale_y_continuous(expand = c(0, 0)) +
  geom_vline(xintercept = 0.58, color = "red")

# SVG
ggsave(plot = p1, filename = "sh_random_knn_5_vs_rho_histogram_12wk.svg", device = "svg", width = 10, height = 10, units = "in", path = "Outputs/knn_analysis/plots/")

# PNG
ggsave(plot = p1, filename = "sh_random_knn_5_vs_rho_histogram_12wk.png", device = "png", width = 10, height = 10, units = "in", path = "Outputs/knn_analysis/plots/")
