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
  purrr, # Function for using map
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

# ----Loading 12 wk data----
twelve_wk_data <- read.csv("Data/Old data/Jacob/12WOcombined_output_MECP2.csv")
twelve_wk_data <- twelve_wk_data %>%
  select(Image, CX..pix., CY..pix., CZ..pix., Condition, Hemishphere, MECP2, Mean) %>%
  rename(
    x = CX..pix.,
    y = CY..pix.,
    z = CZ..pix.,
    mecp2_p = MECP2
  ) %>%
  rename_all(tolower) %>%
  filter(mecp2_p == "P") %>%
  mutate(id = 1:nrow(.)) %>%
  mutate(time = 12) %>%
  mutate(hemishphere = str_extract(image, "(LW|RW)")) %>%
  mutate(hemishphere = ifelse(grepl("LW", hemishphere), gsub("LW", "LH", hemishphere), hemishphere)) %>%
  mutate(hemishphere = ifelse(grepl("RW", hemishphere), gsub("RW", "RH", hemishphere), hemishphere)) %>%
  rename(filename = image)

twelve_wk_data <- twelve_wk_data %>% group_by(filename)
all_imgs <- twelve_wk_data %>% group_split(twelve_wk_data)

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
#' @param x_lims Numeric specifying the x-axis limits of the plots. Default is (0, 0.50).
#' @param plot_mean Should the mean of the random rhos be calculated and then display the difference between it and the observed value? Default is FALSE.
#' @param show_x_value Should the x-intercept value we pass for the obseved value be displayed in red under the vertical line? Default is FALSE.
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
                                  all_search_space = 10000,
                                  calculate_pvalue = TRUE,
                                  p_value_test = "less",
                                  test_type = "Permutation test",
                                  test_direction = "Two-tailed",
                                  x_lims = c(0.0, 0.50),
                                  plot_mean = FALSE,
                                  show_x_value = FALSE) {
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
    # Count how many simulated rhos are as extreme or more extreme than the absolute observed one
    extreme_count <- sum(data$rho >= observed_value)
    message("The number of values as extreme as the observed value of ", observed_value, " are ", extreme_count)
    message("The number of permutations for this calculation is ", num_permutations)
    # Calculate the p-value
    p_value_comp <- extreme_count / all_search_space


    # Print the observed value and p-value to the console
    message(paste0("Observed Value: ", round(observed_value, digits = 4)))
    message(paste0("P-value: ", p_value_comp))

    # Now determining if we need to plot the mean
    if (plot_mean) {
      diff_in_mean <- round(mean(data$rho), digits = 2) - observed_value
      message(paste0("Difference b/t mean and ov: ", diff_in_mean))

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
        ggtitle(paste0("K = ", k_value, " | ", unique(data$condition), " | ", num_permutations, " Permutations | \n", time_point, " WK | ", num_ks, " KNNs Searched", " | Delta Mean: ", diff_in_mean)) +
        xlab(expression(rho)) +
        ylab("Count") +
        scale_y_continuous(expand = c(0, 0)) +
        geom_vline(xintercept = x_intercept, color = line_color) +
        labs(subtitle = paste(test_type, " (", test_direction, ")\np-value: ", p_value_comp, sep = ""))

      # Show x-intercept if need be
      if (show_x_value) {
        p <- p +
          geom_text(aes(x = x_intercept, y = 0, label = round(x_intercept, 4)),
            colour = "red",
            hjust = 1,
            vjust = -1, # Adjust this value for vertical placement of text
            size = 8 * text_scale
          ) # Adjust the size if necessary
      }
    } else {
      # Build the plot with p-value but with no mean plotted
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
        labs(subtitle = paste(test_type, " (", test_direction, ")\np-value: ", p_value_comp, sep = ""))


      if (show_x_value) {
        p <- p +
          geom_text(aes(x = x_intercept, y = 0, label = round(x_intercept, 4)),
            colour = "red",
            hjust = 1,
            vjust = -1, # Adjust this value for vertical placement of text
            size = 8 * text_scale
          ) # Adjust the size if necessary
      }
    }
  } else {
    # Without p-value or mean if not requested
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



    if (show_x_value) {
      p <- p +
        geom_text(aes(x = x_intercept, y = 0, label = round(x_intercept, 4)),
          colour = "red",
          hjust = 1,
          vjust = -1, # Adjust this value for vertical placement of text
          size = 8 * text_scale
        ) # Adjust the size if necessary
    }
  }


  # Calculate the width and height for the plot
  width <- plot_size * dpi
  height <- plot_size * dpi

  # If calculated p-value is less than 2.2e-16 and we are calculating the p-value than add subtitle mentioning raw p-value
  if (calculate_pvalue) {
    if (p_value_comp < 2.2e-16) {
      p <- p + labs(subtitle = paste(test_type, " (", test_direction, ")\np-value: ", round(p_value_comp, digits = 4), sep = ""))
    }
  }

  # Specifying the x-limits for all plots to make comparisons easier
  if (!is.null(x_lims)) {
    p <- p + coord_cartesian(xlim = x_lims)
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
#' Calculate the Euclidean distance between two 3D points
#'
#' This function computes the Euclidean distance between two points in 3D space.
#'
#' @param x1 The x-coordinate of the first point.
#' @param y1 The y-coordinate of the first point.
#' @param z1 The z-coordinate of the first point.
#' @param x2 The x-coordinate of the second point.
#' @param y2 The y-coordinate of the second point.
#' @param z2 The z-coordinate of the second point.
#'
#' @return The Euclidean distance between the two 3D points.
#'
#' @examples
#' euclidean_distance_3d(0, 0, 0, 1, 1, 1)
#' # Output: 1.732051
#'
#' @export
euclidean_distance_3d <- function(x1, y1, z1, x2, y2, z2) {
  sqrt((x2 - x1)^2 + (y2 - y1)^2 + (z2 - z1)^2)
}

# ----Function to perform each semi-random KNN analysis----
#' Process a single iteration of the semi-random KNN analysis
#'
#' @param iteration The current iteration number.
#' @param max_iterations The maximum number of iterations to run (default is 1000).
#' @param time_point_to_use The time point to use (default is 6).
#' @param all_imgs A list of data frames containing image data.
#' @param num_nearest The number of nearest neighbors to consider.
#' @param k_neighbors The number of neighbors to sample from the nearest neighbors.
#' @param apply_intensity Logical, indicating whether to apply intensity information (default is TRUE).
#' @param p_only Logical, indicating whether to include only "P" labeled neighbors (default is TRUE).
#' @param calc_mean_intensity Logical, indicating whether to calculate the mean intensity (default is TRUE).
#' @param conditions A character vector of conditions to consider (default is c("NW", "SW", "NH", "SH")).
#'
#' @return A list of correlation values for each condition.
#'
#' @export

process_iteration <- function(iteration, max_iterations = 1000, 
                              time_point_to_use = 6, 
                              all_imgs, num_nearest, 
                              k_neighbors = 5,
                              apply_intensity = TRUE,
                              p_only = TRUE, 
                              calc_mean_intensity = TRUE, 
                              conditions = c("NW", "SW", "NH", "SH")) {
  random_seed <- runif(1, max_iterations, 100000)
  write.table(random_seed,  
               file= paste0("Outputs/knn_analysis/random_seeds_for_semi_random_analysis/random_seeds_for_",time_point_to_use,"_wk_semi-random.csv"), 
               append = T, 
               sep=',', 
               row.names=F, 
               col.names=F )
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
        k_nearest_points <- num_nearest_points[sampled_indices, ]
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
      correlation_values[[cond]] <- current_df
    }
  }
  return(correlation_values)
}


# ----Function to perform purely random KNN analysis per image----
process_single_image_iteration <- function(iteration,
                                           max_iterations = 1000, 
                                           time_point_to_use = 6, 
                                           img = NULL, 
                                           k_neighbors = 5,
                                           apply_intensity = TRUE,
                                           p_only = TRUE, 
                                           calc_mean_intensity = TRUE, 
                                           conditions = NULL) {
  
  random_seed <- runif(1, max_iterations, 100000)
  
  set.seed(random_seed)
  
  iteration_data <- list()
  counter <- 1
  
  df <- img
  df$knn_5_label <- rep(0, nrow(df))
  
  for (r in 1:nrow(df)) {
    observation <- df[r, ]
    other_points <- df[-r, ]
    sampled_indices <- sample(1:nrow(other_points), k_neighbors, replace = TRUE)
    k_nearest_points <- other_points[sampled_indices, ]
    
    intensity_info <- k_nearest_points
    add_p_score <- sum((intensity_info$mecp2_p == "P")) / k_neighbors
    observation$knn_5_label <- add_p_score
    neighbor_means <- c()
    
    if (apply_intensity) {
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
  
  combined_data <- dplyr::bind_rows(iteration_data)
  return(combined_data)
}


# ----6 Week Semi-random parallel code----
# Set the number of CPU cores to use
num_cores <- 11

# Initialize parallel backend
cl <- makeCluster(num_cores)
registerDoParallel(cl)


# With progress bar
total_iterations <- 1
all_results <- list()
counter <- 1
search_space <- seq(5, 50, 5)

# Define the total number of tasks
# total_tasks <- length(search_space) * 1000

total_tasks <- length(search_space) * 10

# Create a progress bar
pb <- txtProgressBar(min = 0, max = total_tasks, style = 3)

for (s in search_space) {
  for (i in 1:10) {
    results <- foreach(iteration = 1:total_iterations, .combine = c) %dopar% {
      process_iteration(iteration, max_iterations = total_iterations, time_point_to_use = 6, all_imgs = all_imgs, num_nearest = s, k_neighbors = 5, apply_intensity = TRUE, p_only = TRUE, calc_mean_intensity = TRUE, conditions = TRUE)
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

# ----Binding the 6 wk results together----
# This list will store the combined dataframes
all_combined_dfs <- list()

# Loop over the six_wk_sr_res list
for (name in names(six_wk_sr_res)) {
  # For each condition (NW/NH/SW/SH), bind rows and add the condition
  for (condition in names(six_wk_sr_res[[name]])) {
    df <- six_wk_sr_res[[name]][[condition]][1,]
    df$Source <- name
    df$condition <- condition
    all_combined_dfs[[paste(name, condition, sep = "_")]] <- df
  }
}

# Now bind all dataframes into a single one
final_df <- bind_rows(all_combined_dfs)


# Extract the iteration numbers from the Source column
iteration_numbers <- as.numeric(str_extract(final_df$Source, "(?<=Iteration )[0-9]+"))

# Identify the largest iteration number
max_iteration <- max(iteration_numbers, na.rm = TRUE)

# Then divide this number by the length of search_space to get number of iterations performed per search value
max_iteration <- max_iteration / length(search_space)

# Update the num_iterations column with the largest iteration number
final_df$num_iterations <- max_iteration

# ----Plotting 6 wk semi-random 1000 rep results----
conds <- c("NW", "NH", "SW", "SH")
x_ints <- c(0.27, 0.44, 0.52, 0.58)


all_plots <- list()
current_combo_plot <- list()
counter <- 1
for (cond in conds) {
  for (x in x_ints[counter]) {
    for (s in search_space) {
      current_df <- filter(final_df, condition == cond)
      current_df <- filter(current_df, search_space == s)
      current_plot <- plot_permutation_test(
        data = current_df,
        filename = paste0(tolower(cond), "_", s, "_nn_6wk"),
        x_intercept = x,
        calculate_pvalue = TRUE,
        plot_size = 2.5,
        base_plot_size = 8,
        dpi = 600,
        x_lims = c(0, 0.60)
      )
      all_plots[[paste0(cond, "_", x, "_", s)]] <- current_plot
      current_combo_plot[[paste0(cond, "_", x, "_", s)]] <- current_plot
    }
  }
  current_big_plot <- ggarrange(plotlist = current_combo_plot, ncol = 2, nrow = 5, labels = "AUTO", align = "hv", font.label = list(size = 16, face = "bold"))
  ggsave(plot = current_big_plot, device = "png", path = "Outputs/knn_analysis/plots/", filename = paste0(tolower(cond), "_all_serch_space_combo_plot_6wk_new_p_value_code.png"), width = 8, height = 8, units = "in", dpi = 600, bg = "white")
  current_big_plot <- list()
  current_combo_plot <- list()
  counter <- counter + 1
}


# ----12 Week Semi-random Parallel Code----
# Set the number of CPU cores to use (you can adjust this based on your system)
num_cores <- 11

# Initialize parallel backend
cl <- makeCluster(num_cores)
registerDoParallel(cl)


# With progress bar
total_iterations <- 1
all_results <- list()
counter <- 1
search_space <- seq(5, 20, 5)

# Define the total number of tasks
total_tasks <- length(search_space) * 1000

# Create a progress bar
pb <- txtProgressBar(min = 0, max = total_tasks, style = 3)

for (s in search_space) {
  for (i in 1:1000) {
    results <- foreach(iteration = 1:total_iterations, .combine = c) %dopar% {
      process_iteration(iteration, max_iterations = total_iterations, all_imgs = all_imgs, num_nearest = s, k_neighbors = 5, time_point_to_use = 12, apply_intensity = TRUE, p_only = TRUE, calc_mean_intensity = TRUE, conditions = TRUE)
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


# ----Binding the 12 wk results together----
# This list will store the combined dataframes
all_combined_dfs <- list()

# Loop over the all_results list
for (name in names(all_results)) {
  # For each condition (NW/NH/SW/SH), bind rows and add the condition
  for (condition in names(all_results[[name]])) {
    df <- all_results[[name]][[condition]][1,]
    df$Source <- name
    df$condition <- condition
    all_combined_dfs[[paste(name, condition, sep = "_")]] <- df
  }
}

# Now bind all dataframes into a single one
final_df <- bind_rows(all_combined_dfs)


# Extract the iteration numbers from the Source column
iteration_numbers <- as.numeric(str_extract(final_df$Source, "(?<=Iteration )[0-9]+"))

# Identify the largest iteration number
max_iteration <- max(iteration_numbers, na.rm = TRUE)

# Then divide this number by the length of search_space to get number of iterations performed per search value
max_iteration <- max_iteration / length(search_space)

# Update the num_iterations column with the largest iteration number
final_df$num_iterations <- max_iteration

# ----Plotting 12 wk semi-random 1000 rep results----
conds <- c("NW", "NH", "SW", "SH")
x_ints <- c(0.23, 0.17, 0.34, 0.26)


all_plots <- list()
current_combo_plot <- list()
counter <- 1
for (cond in conds) {
  for (x in x_ints[counter]) {
    for (s in search_space) {
      current_df <- filter(final_df, condition == cond)
      current_df <- filter(current_df, search_space == s)
      current_plot <- plot_permutation_test(
        data = current_df,
        filename = paste0(tolower(cond), "_", s, "_nn_12_wk"),
        x_intercept = x,
        calculate_pvalue = TRUE,
        plot_size = 2.5,
        base_plot_size = 8,
        dpi = 600,
        x_lims = c(0, 0.60),
        plot_mean = TRUE
      )
      all_plots[[paste0(cond, "_", x, "_", s)]] <- current_plot
      current_combo_plot[[paste0(cond, "_", x, "_", s)]] <- current_plot
    }
  }
  current_big_plot <- ggarrange(plotlist = current_combo_plot, ncol = 2, nrow = 2, labels = "AUTO", align = "hv", font.label = list(size = 16, face = "bold"))
  ggsave(plot = current_big_plot, device = "png", path = "Outputs/knn_analysis/plots/", filename = paste0(tolower(cond), "_all_serch_space_combo_plot_12_wk.png"), width = 8, height = 8, units = "in", dpi = 600, bg = "white")
  current_big_plot <- list()
  current_combo_plot <- list()
  counter <- counter + 1
}


# ----Completely random 6wk code by both image and condition----
finished_imgs <- list()
all_iterations <- list()
counter <- 1

# Create a progress bar
pb <- txtProgressBar(min = 0, max = 1000*length(all_imgs), style = 3)

for (i in 1:length(all_imgs)) {
  current_img <- all_imgs[[i]]
  for (i in 1:1000) {
    current_data <- process_single_image_iteration(img = current_img)
    finished_imgs[[paste0("Iteration_",counter)]] <- current_data
    
    # Update the progress bar
    setTxtProgressBar(pb, counter)
    
    counter <- counter + 1
  }
}


# Close the progress bar
close(pb)

# Now calculating the rhos for each of the images
# Split the list into a list of data frames
list_of_dfs <- map(seq(1, length(finished_imgs), by = 1000), function(start_index) {
  end_index <- start_index + 999
  tibble(elements = finished_imgs[start_index:end_index])
})


# Function to compute Spearman correlation for a given iteration of an image
calculate_correlation <- function(iteration_data) {
  vec1 <- as.numeric(iteration_data$mean)
  vec2 <- as.numeric(iteration_data$positive_neighborhood_mean)
  
  # Check if vec1 and vec2 have data
  if(length(vec1) == 0 || length(vec2) == 0) {
    return(NA)
  }
  
  cor_val <- cor.test(vec1, vec2, method = "spearman", use = "complete.obs")$estimate
  return_df <- data.frame(condition = iteration_data$condition[1],
                          rho = as.vector(unname(cor_val)), 
                          time_point = iteration_data$time_point[1], 
                          k_value = iteration_data$k_value[1], 
                          search_space = iteration_data$search_space[1], 
                          num_iterations = iteration_data$num_iterations[1])
  return(return_df)
}


# Process each image
results <- map(seq_along(list_of_dfs), function(image_num) {
  
  image_df <- list_of_dfs[[image_num]]
  
  # Print the current image number to verify
  print(paste0("Processing image ", image_num))
  
  # Process each iteration for the current image
  iterations_results <- map(1:1000, function(iter_num) {
    
    # Extract the data for the current iteration
    iteration_data <- image_df[[1]][[iter_num]]
    iteration_data$time_point <- rep(6, nrow(iteration_data))
    iteration_data$k_value <- rep(5, nrow(iteration_data))
    iteration_data$search_space <- "max_randomness"
    iteration_data$num_iterations <- 1000
    
    
    # Calculate the Spearman correlation for the current iteration
    cor_val <- calculate_correlation(iteration_data)
    return(cor_val)
  })
  
  return(iterations_results)
})



# Now binding all the iterations together per image
bound_results <- map(results, bind_rows)


# Now plotting all of the 6 wk results
# conds <- c("NW", "NH", "SW", "SH")
# x_ints <- c(0.27, 0.44, 0.52, 0.58)

condition_to_x_intercept <- function(condition) {
  lookup <- list(
    "NW" = 0.27,
    "NH" = 0.44,
    "SW" = 0.52,
    "SH" = 0.58
  )
  
  return(lookup[[condition]])
}


# Generate filenames based on the image number
filenames <- paste0("Image_", 1:length(bound_results))

# Iterate over each dataframe and call the plot_permutation_test function
plots <- map2(bound_results, filenames, function(data, fname) {
  # Extract the unique condition from the current dataframe
  current_condition <- unique(data$condition)
  
  # Get the corresponding x_intercept for the current condition
  x_intercept_value <- condition_to_x_intercept(current_condition)
  
  # Call the plot_permutation_test function with the determined x_intercept
  plot_permutation_test(data = data, filename = fname, x_intercept = x_intercept_value,
                        all_search_space = 1000, 
                        calculate_pvalue = TRUE, 
                        plot_size = 2.5,
                        base_plot_size = 8, 
                        dpi = 600, 
                        x_lims = c(-0.60,0.60),
                        plot_mean = TRUE, 
                        show_x_value = FALSE)
})



# Combo plot
condition_list <- sapply(plots, function(p) {
  unique(p$data$condition)
})

order_indices <- order(factor(condition_list, levels = c("NW", "NH", "SW", "SH")))
ordered_plots <- plots[order_indices]

combo_plot <- ggarrange(plotlist = ordered_plots, ncol = 8, nrow = 7)
ggsave(combo_plot, bg = "white", path = "Outputs/knn_analysis/plots/", filename = "combo_plot_6wk_by_img_total_random.png", device = "png", width = 20, height = 20, units = "in", dpi = 600)
ggsave(combo_plot, bg = "white", path = "Outputs/knn_analysis/plots/", filename = "combo_plot_6wk_by_img_total_random.svg", device = "svg", width = 20, height = 20, units = "in", dpi = 600)

# ----Completely random 12wk code by both image and condition----
finished_imgs <- list()
all_iterations <- list()
counter <- 1

# Create a progress bar
pb <- txtProgressBar(min = 0, max = 1000*length(all_imgs), style = 3)

for (i in 1:length(all_imgs)) {
  current_img <- all_imgs[[i]]
  for (i in 1:1000) {
    current_data <- process_single_image_iteration(img = current_img, time_point_to_use = 12)
    finished_imgs[[paste0("Iteration_",counter)]] <- current_data
    
    # Update the progress bar
    setTxtProgressBar(pb, counter)
    
    counter <- counter + 1
  }
}


# Close the progress bar
close(pb)

# Now calculating the rhos for each of the images
# Split the list into a list of data frames
list_of_dfs <- map(seq(1, length(finished_imgs), by = 1000), function(start_index) {
  end_index <- start_index + 999
  tibble(elements = finished_imgs[start_index:end_index])
})


# Function to compute Spearman correlation for a given iteration of an image
calculate_correlation <- function(iteration_data) {
  vec1 <- as.numeric(iteration_data$mean)
  vec2 <- as.numeric(iteration_data$positive_neighborhood_mean)
  
  # Check if vec1 and vec2 have data
  if(length(vec1) == 0 || length(vec2) == 0) {
    return(NA)
  }
  
  cor_val <- cor.test(vec1, vec2, method = "spearman", use = "complete.obs")$estimate
  return_df <- data.frame(condition = iteration_data$condition[1],
                          rho = as.vector(unname(cor_val)), 
                          time_point = iteration_data$time_point[1], 
                          k_value = iteration_data$k_value[1], 
                          search_space = iteration_data$search_space[1], 
                          num_iterations = iteration_data$num_iterations[1])
  return(return_df)
}


# Process each image
results <- map(seq_along(list_of_dfs), function(image_num) {
  
  image_df <- list_of_dfs[[image_num]]
  
  # Print the current image number to verify
  print(paste0("Processing image ", image_num))
  
  # Process each iteration for the current image
  iterations_results <- map(1:1000, function(iter_num) {
    
    # Extract the data for the current iteration
    iteration_data <- image_df[[1]][[iter_num]]
    iteration_data$time_point <- rep(12, nrow(iteration_data))
    iteration_data$k_value <- rep(5, nrow(iteration_data))
    iteration_data$search_space <- "max_randomness"
    iteration_data$num_iterations <- 1000
    
    
    # Calculate the Spearman correlation for the current iteration
    cor_val <- calculate_correlation(iteration_data)
    return(cor_val)
  })
  
  return(iterations_results)
})



# Now binding all the iterations together per image
bound_results <- map(results, bind_rows)


# Now plotting all of the 12 wk results
# conds <- c("NW", "NH", "SW", "SH")
# x_ints <- c(0.23, 0.17, 0.34, 0.26)

condition_to_x_intercept <- function(condition) {
  lookup <- list(
    "NW" = 0.23,
    "NH" = 0.17,
    "SW" = 0.35,
    "SH" = 0.26
  )
  
  return(lookup[[condition]])
}


# Generate filenames based on the image number
filenames <- paste0("Image_", 1:length(bound_results))

# Iterate over each dataframe and call the plot_permutation_test function
plots <- map2(bound_results, filenames, function(data, fname) {
  # Extract the unique condition from the current dataframe
  current_condition <- unique(data$condition)
  
  # Get the corresponding x_intercept for the current condition
  x_intercept_value <- condition_to_x_intercept(current_condition)
  
  # Call the plot_permutation_test function with the determined x_intercept
  plot_permutation_test(data = data, filename = fname, x_intercept = x_intercept_value,
                        all_search_space = 1000, 
                        calculate_pvalue = TRUE, 
                        plot_size = 2.5,
                        base_plot_size = 8, 
                        dpi = 600, 
                        x_lims = c(-0.60,0.60),
                        plot_mean = TRUE, 
                        show_x_value = FALSE)
})



# Combo plot
condition_list <- sapply(plots, function(p) {
  unique(p$data$condition)
})

order_indices <- order(factor(condition_list, levels = c("NW", "NH", "SW", "SH")))
ordered_plots <- plots[order_indices]

combo_plot <- ggarrange(plotlist = ordered_plots, ncol = 5, nrow = 10)
ggsave(combo_plot, bg = "white", path = "Outputs/knn_analysis/plots/", filename = "combo_plot_12wk_by_img_total_random.png", device = "png", width = 20, height = 20, units = "in", dpi = 600)
ggsave(combo_plot, bg = "white", path = "Outputs/knn_analysis/plots/", filename = "combo_plot_12wk_by_img_total_random.svg", device = "svg", width = 20, height = 20, units = "in", dpi = 600)

# ----Completely random 6 wk KNN analysis by just condition [like initial analysis]----
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
    
    # Instead of calculating distances, sample random indices. Initially replace was set to FALSE 
    sampled_indices <- sample(1:nrow(other_points), k_neighbors, replace = TRUE)
    k_nearest_points <- other_points[sampled_indices, ]
    intensity_info <- other_points[sampled_indices, ]
    
    add_p_score <- sum((intensity_info$mecp2_p == "P")) / k_neighbors
    observation$knn_5_label <- add_p_score
    all_p_scores[[counter]] <- add_p_score
    
    # Calculate distances with intensity information
    if (apply_intensity) {
      # No need to filter by k_nearest_indices since we're directly sampling them
      distances_updated <- rescale(intensity_info$mean, to = c(0, 1))
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
        mean_of_neighbors <- mean(intensity_info$mean)
        neighbor_means <- c(neighbor_means, mean_of_neighbors)
        observation$positive_neighborhood_mean <- mean_of_neighbors
        all_data[[counter]] <- observation
        counter <- counter + 1
      }
    }
  }
}
test_df <- bind_rows(all_data)


# Total number of iterations
num_iterations <- 1000

# Initialize the progress bar
pb <- txtProgressBar(min = 0, max = num_iterations, style = 3)

all_iterations_results <- list()

for (iter in 1:num_iterations) {
  
  all_p_scores <- list()
  all_data <- list()
  counter <- 1
  
  for (i in 1:length(all_imgs)) {
    df <- all_imgs[[i]]
    df$knn_5_label <- rep(0, nrow(df))
    
    for (r in 1:nrow(df)) {
      observation <- df[r, ]
      other_points <- df[-r, ]
      
      # Sample random indices
      sampled_indices <- sample(1:nrow(other_points), k_neighbors)
      k_nearest_points <- other_points[sampled_indices, ]
      intensity_info <- other_points[sampled_indices, ]
      
      add_p_score <- sum((intensity_info$mecp2_p == "P")) / k_neighbors
      observation$knn_5_label <- add_p_score
      all_p_scores[[counter]] <- add_p_score
      
      # Calculate distances with intensity information
      if (apply_intensity) {
        distances_updated <- rescale(intensity_info$mean, to = c(0, 1))
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
          mean_of_neighbors <- mean(intensity_info$mean)
          observation$positive_neighborhood_mean <- mean_of_neighbors
          all_data[[counter]] <- observation
          counter <- counter + 1
        }
      }
    }
  }
  
  all_iterations_results[[iter]] <- all_data
  
  # Update the progress bar
  setTxtProgressBar(pb, iter)
}


# Close the progress bar after completion
close(pb)

all_cond_data <- bind_rows(all_iterations_results)


# Now calculating each of the rhos
conds <- unique(all_cond_data$condition)
num_iterations <- 1000
all_results <- list()


for (cond in conds) {
  current_cond <- filter(all_cond_data, condition == cond)
  num_records <- nrow(current_cond)
  records_per_iteration <- num_records / num_iterations
  
  cor_results <- numeric(num_iterations)  # To store the correlation coefficients
  
  for (i in 1:num_iterations) {
    start_idx <- (i - 1) * records_per_iteration + 1
    # message(paste0("Start id for cond: ", cond, " iteration ", i, " is: ", start_idx))
    end_idx <- i * records_per_iteration
    # message(paste0("End id for cond: ", cond, "iteration ", i, "is: ", end_idx))
    subset_df <- current_cond[start_idx:end_idx, ]
    
    # Check if there are enough records in both vectors to compute correlation
    if(length(unique(subset_df$positive_neighborhood_mean)) > 1 && length(unique(subset_df$mean)) > 1) {
      cor_val <- cor.test(subset_df$mean, subset_df$positive_neighborhood_mean, method = "spearman", use = "complete.obs")$estimate
      cor_results[i] <- cor_val
    } else {
      cor_results[i] <- NA
    }
  }
  
  all_results[[paste0(cond)]] <- data.frame(condition = rep(cond, length(cor_results)), 
                                            rho = cor_results,
                                            time_point = rep(6, length(cor_results)),
                                            k_value = rep(5, length(cor_results)),
                                            num_iterations = rep(num_iterations, length(cor_results)),
                                            search_space = rep("max_randomness", length(cor_results)))
  
}

all_plots <- list()
conds <- c("NW", "NH", "SW", "SH")
x_ints <- c(0.27, 0.44, 0.52, 0.58)

# Re-order list if needed
all_results<- all_results[c("NW", "NH", "SW", "SH")]

for (cond in 1:length(conds)) {
  current_plot <- plot_permutation_test(data = all_results[[cond]], 
                                        filename = paste0(tolower(conds[cond]),"_1000_random_perm_k_5_6wk_with_replacement"), 
                                        plot_size = 2.5, 
                                        calculate_pvalue = TRUE, 
                                        plot_mean = TRUE, 
                                        x_intercept = x_ints[cond],
                                        x_lims = c(-0.60, 0.60))
  
  all_plots[[paste0(cond)]] <- current_plot
}

combo_plot <- ggarrange(plotlist = all_plots,
                        ncol = 2,
                        nrow = 2, 
                        labels = "AUTO", 
                        font.label = list(size = 18, face = "bold"),
                        align = "hv")

ggsave(combo_plot, device = "png", path = "Outputs/knn_analysis/plots/", filename = "combo_plot_all_conditions_6wk_data_with_replacement.png", units = "in", width = 8, height = 8, bg = "white")



# ----Completely random 12 wk KNN analysis by just condition [like initial analysis]----
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
    
    # Instead of calculating distances, sample random indices
    sampled_indices <- sample(1:nrow(other_points), k_neighbors, replace = TRUE)
    k_nearest_points <- other_points[sampled_indices, ]
    intensity_info <- other_points[sampled_indices, ]
    
    add_p_score <- sum((intensity_info$mecp2_p == "P")) / k_neighbors
    observation$knn_5_label <- add_p_score
    all_p_scores[[counter]] <- add_p_score
    
    # Calculate distances with intensity information
    if (apply_intensity) {
      # No need to filter by k_nearest_indices since we're directly sampling them
      distances_updated <- rescale(intensity_info$mean, to = c(0, 1))
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
        mean_of_neighbors <- mean(intensity_info$mean)
        neighbor_means <- c(neighbor_means, mean_of_neighbors)
        observation$positive_neighborhood_mean <- mean_of_neighbors
        all_data[[counter]] <- observation
        counter <- counter + 1
      }
    }
  }
}
test_df <- bind_rows(all_data)


# Total number of iterations
num_iterations <- 1000

# Initialize the progress bar
pb <- txtProgressBar(min = 0, max = num_iterations, style = 3)

all_iterations_results <- list()

for (iter in 1:num_iterations) {
  
  all_p_scores <- list()
  all_data <- list()
  counter <- 1
  
  for (i in 1:length(all_imgs)) {
    df <- all_imgs[[i]]
    df$knn_5_label <- rep(0, nrow(df))
    
    for (r in 1:nrow(df)) {
      observation <- df[r, ]
      other_points <- df[-r, ]
      
      # Sample random indices
      sampled_indices <- sample(1:nrow(other_points), k_neighbors, replace = TRUE)
      k_nearest_points <- other_points[sampled_indices, ]
      intensity_info <- other_points[sampled_indices, ]
      
      add_p_score <- sum((intensity_info$mecp2_p == "P")) / k_neighbors
      observation$knn_5_label <- add_p_score
      all_p_scores[[counter]] <- add_p_score
      
      # Calculate distances with intensity information
      if (apply_intensity) {
        distances_updated <- rescale(intensity_info$mean, to = c(0, 1))
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
          mean_of_neighbors <- mean(intensity_info$mean)
          observation$positive_neighborhood_mean <- mean_of_neighbors
          all_data[[counter]] <- observation
          counter <- counter + 1
        }
      }
    }
  }
  
  all_iterations_results[[iter]] <- all_data
  
  # Update the progress bar
  setTxtProgressBar(pb, iter)
}


# Close the progress bar after completion
close(pb)

all_cond_data <- bind_rows(all_iterations_results)


# Now calculating each of the rhos
conds <- unique(all_cond_data$condition)
num_iterations <- 1000
all_results <- list()


for (cond in conds) {
  current_cond <- filter(all_cond_data, condition == cond)
  num_records <- nrow(current_cond)
  records_per_iteration <- num_records / num_iterations
  
  cor_results <- numeric(num_iterations)  # To store the correlation coefficients
  
  for (i in 1:num_iterations) {
    start_idx <- (i - 1) * records_per_iteration + 1
    end_idx <- i * records_per_iteration
    subset_df <- current_cond[start_idx:end_idx, ]
    
    # Check if there are enough records in both vectors to compute correlation
    if(length(unique(subset_df$positive_neighborhood_mean)) > 1 && length(unique(subset_df$mean)) > 1) {
      cor_val <- cor.test(subset_df$mean, subset_df$positive_neighborhood_mean, method = "spearman", use = "complete.obs")$estimate
      cor_results[i] <- cor_val
    } else {
      cor_results[i] <- NA
    }
  }
  
  all_results[[paste0(cond)]] <- data.frame(condition = rep(cond, length(cor_results)), 
                                            rho = cor_results,
                                            time_point = rep(12, length(cor_results)),
                                            k_value = rep(5, length(cor_results)),
                                            num_iterations = rep(num_iterations, length(cor_results)),
                                            search_space = rep("max_randomness", length(cor_results)))
  
}

all_plots <- list()
conds <- c("NW", "NH", "SW", "SH")
x_ints <- c(0.23, 0.17, 0.34, 0.26)

# Re-order list if needed
all_results<- all_results[c("NW", "NH", "SW", "SH")]

for (cond in 1:length(conds)) {
  current_plot <- plot_permutation_test(data = all_results[[cond]], 
                                        filename = paste0(tolower(conds[cond]),"_1000_random_perm_k_5_12wk"), 
                                        plot_size = 2.5, 
                                        calculate_pvalue = TRUE, 
                                        plot_mean = TRUE, 
                                        x_intercept = x_ints[cond],
                                        x_lims = c(0.0, 0.35))
  
  all_plots[[paste0(cond)]] <- current_plot
}

combo_plot <- ggarrange(plotlist = all_plots,
                        ncol = 2,
                        nrow = 2, 
                        labels = "AUTO", 
                        font.label = list(size = 18, face = "bold"),
                        align = "hv")

ggsave(combo_plot, device = "png", path = "Outputs/knn_analysis/plots/", filename = "combo_plot_all_conditions_12wk_data.png", units = "in", width = 8, height = 8, bg = "white")
