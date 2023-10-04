# ----Loading packages----
pacman::p_load(
  car, # Tools for performing specific statistical tests
  caret, # Tools for training and evaluating predictive models
  class, # Functions for nearest neighbor classification
  deldir, # Delaunay triangulation and Dirichlet (Voronoi) tesselation
  doParallel, # Additional parallel functions
  effectsize, # Functions for calculating standardized effect sizes
  foreach, # Functions for parallel computing
  ggbreak, # Function for allowing creating of breaks in plots
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
  # filter(mecp2_p == "P") %>%
  mutate(id = 1:nrow(.))

six_wk_data <- six_wk_data %>% group_by(filename)
all_imgs <- six_wk_data %>% group_split(six_wk_data)


# ----Loading revised 6 wk data----
six_wk_data <- read.csv("Data/Old data/Jacob/3D_MECP2_6WO_combined_dataset_revised.csv")
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
  # Comment out this line if doing n & p cell analysis
  # filter(mecp2_p == "P") %>%
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
  # Comment out this line if doing n & p cell analysis
  # filter(mecp2_p == "P") %>%
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
                                  all_search_space = 1000,
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
p_only <- FALSE
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

# Save the 6 wk dataframe used to calculate correlation coefs
saveRDS(test_df, paste0("Outputs/knn_analysis/data/six_wk_df_to_calculate_correlation_rhos_revised_with_only_p_mecp2_samples_z_corrected_data_", Sys.time(), ".rds"))

# Total number of iterations
num_iterations <- 1000

# Initialize the progress bar
pb <- progress_bar$new(total = num_iterations, format = "[:bar] :percent :current/:total :eta")

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
  pb$tick()
}


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
# x_ints <- c(0.27, 0.44, 0.52, 0.58)
x_ints <- c(0.37, 0.43, 0.50, 0.51)

# Re-order list if needed
all_results<- all_results[c("NW", "NH", "SW", "SH")]

for (cond in 1:length(conds)) {
  current_plot <- plot_permutation_test(data = all_results[[cond]], 
                                        filename = paste0(tolower(conds[cond]),"_1000_random_perm_k_5_6wk_with_replacement_corrected_p_value_revised_data_p_and_n_neighbors_included_from_initial_data_finished"), 
                                        plot_size = 2.5, 
                                        calculate_pvalue = TRUE, 
                                        plot_mean = TRUE, 
                                        all_search_space = 1000,
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

ggsave(combo_plot, device = "png", path = "Outputs/knn_analysis/plots/", filename = "combo_plot_all_conditions_6wk_data_with_replacement_corrected_p_value_revised_p_and_n_neighbors_included_from_initial_data_finished.png", units = "in", width = 8, height = 8, bg = "white")



# ----Completely random 12 wk KNN analysis by just condition [like initial analysis]----
all_data <- list()
all_p_scores <- list()
counter <- 1
apply_intensity <- TRUE
p_only <- FALSE
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
pb <- progress_bar$new(total = num_iterations, format = "[:bar] :percent :current/:total :eta")


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
  pb$tick()
}


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
x_ints <- c(0.32, 0.20, 0.24, 0.22)
# x_ints <- c(0.23, 0.20, 0.34, 0.28)

# Re-order list if needed
all_results<- all_results[c("NW", "NH", "SW", "SH")]

for (cond in 1:length(conds)) {
  current_plot <- plot_permutation_test(data = all_results[[cond]], 
                                        filename = paste0(tolower(conds[cond]),"_1000_random_perm_k_5_12wk_corrected_p_value_n_and_p_neighbors_only_included_from_initial_data_with_updated_observed_rho"), 
                                        plot_size = 2.5, 
                                        calculate_pvalue = TRUE, 
                                        plot_mean = TRUE, 
                                        all_search_space = 1000,
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

ggsave(combo_plot, device = "png", path = "Outputs/knn_analysis/plots/", filename = "combo_plot_all_conditions_12wk_data_corrected_p_value_n_and_p_neighbors_only_included_from_initial_data_finished.png", units = "in", width = 8, height = 8, bg = "white")


# ----True K=5 neighbors 6 wk data----
all_data <- list()
all_p_scores <- list()
counter <- 1
apply_intensity <- TRUE
p_only <- FALSE
calc_mean_intensity <- TRUE
neighbor_means <- c()


k_neighbors <- 5
for (i in 1:length(all_imgs)) {
  df <- all_imgs[[i]]
  df$knn_5_label <- rep(0, nrow(df))
  
  for (r in 1:nrow(df)) {
    observation <- df[r, ]
    other_points <- df[-r, ]
    
    # Calculating distances
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

# Save the 6 wk dataframe used to calculate correlation coefs
saveRDS(test_df, paste0("Outputs/knn_analysis/data/six_wk_df_to_calculate_correlation_rhos_revised_with_both_n_and_p_mecp2_samples_true_k_5_neighbors", Sys.time(), ".rds"))

# Total number of iterations
num_iterations <- 1000

# Initialize the progress bar
pb <- progress_bar$new(total = num_iterations, format = "[:bar] :percent :current/:total :eta")

all_iterations_results <- list()

for (iter in 1:num_iterations) {
  
  all_p_scores <- list()
  all_data <- list()
  counter <- 1
  k <- 5
  
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
  pb$tick()
}




all_cond_data <- bind_rows(all_iterations_results)


# ----Remaking scatter plots with the n and p neighbors in the data included for 6 wk----
plot_type <- "png"
time_point <- "6wk"
fit_line <- "reg.line"
stat_used <- "spearman"
plots <- list()
k_n <- 5
counter <- 1


mean_neighbor_df_simplified <- test_df
mean_neighbor_df_simplified <- mean_neighbor_df_simplified %>%
  mutate(filename = match(filename, unique(filename))) %>%
  rename(Image = filename)


for (c in unique(mean_neighbor_df_simplified$condition)) {
  current_cond <- filter(mean_neighbor_df_simplified, condition == c)
  current_cond$Image <- as.factor(current_cond$Image)
  # current_cond <- filter(current_cond, mecp2_p == "P")
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
    labs(color = "Image", fill = "Image")+
    coord_cartesian(xlim = c(0, 8000), ylim = c(0, 3000))
  
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
       device = plot_type, path = "Outputs/knn_analysis/plots/",
       filename = paste0("scatter_plot_3d_combo_p_cells_included_", "_", time_point, ".", plot_type),
       width = 14, height = 14, units = "in", dpi = 600, bg = "white"
)


# Saving each individual plot
for (i in 1:length(plots)) {
  ggsave(
    plot = plots[[i]],
    filename = paste0("scatter_plot_3d_p_cells_included_", i, "_", time_point, ".", plot_type),
    width = 10, height = 10, units = "in", dpi = 600, bg = "white", device = plot_type, path = "Outputs/knn_analysis/plots/"
  )
}

# ----True K=5 neighbors 12 wk data----
all_data <- list()
all_p_scores <- list()
counter <- 1
apply_intensity <- TRUE
p_only <- FALSE
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

# Save the 12 wk dataframe used to calculate correlation coefs
saveRDS(test_df, paste0("Outputs/knn_analysis/data/twelve_wk_df_to_calculate_correlation_rhos_revised_with_both_n_and_p_mecp2_samples_true_k_5_neighbors", Sys.time(), ".rds"))

# Total number of iterations
num_iterations <- 1000

# Initialize the progress bar
pb <- progress_bar$new(total = num_iterations, format = "[:bar] :percent :current/:total :eta")

all_iterations_results <- list()

for (iter in 1:num_iterations) {
  
  all_p_scores <- list()
  all_data <- list()
  counter <- 1
  k <- 5
  
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
  pb$tick()
}

all_cond_data <- bind_rows(all_iterations_results)




# ----Remaking scatter plots with the n and p neighbors in the data included for 12 wk----
plot_type <- "png"
time_point <- "12wk"
fit_line <- "reg.line"
stat_used <- "spearman"
plots <- list()
k_n <- 5
counter <- 1


mean_neighbor_df_simplified <- test_df
mean_neighbor_df_simplified <- mean_neighbor_df_simplified %>%
  mutate(filename = match(filename, unique(filename)))%>%
  rename(Image = filename)



for (c in unique(mean_neighbor_df_simplified$condition)) {
  current_cond <- filter(mean_neighbor_df_simplified, condition == c)
  current_cond$Image <- as.factor(current_cond$Image)
  # current_cond <- filter(current_cond, mecp2_p == "P")
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
    ylab = "Mean Neighborhood Intensity",
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
    labs(color = "Image", fill = "Image")+
    coord_cartesian(xlim = c(0, 8000), ylim = c(0,3000))
  
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
       device = plot_type, path = "Outputs/knn_analysis/plots/",
       filename = paste0("scatter_plot_3d_combo_p_and_n_cells_included_", "_", time_point, ".", plot_type),
       width = 14, height = 14, units = "in", dpi = 600, bg = "white"
)


# Saving each individual plot
for (i in 1:length(plots)) {
  ggsave(
    plot = plots[[i]],
    filename = paste0("scatter_plot_3d_p_and_n_cells_included_", i, "_", time_point, ".", plot_type),
    width = 10, height = 10, units = "in", dpi = 600, bg = "white", device = plot_type, path = "Outputs/knn_analysis/plots/"
  )
}

# ----True K neighbors for 12 wk data for n and p cells together----
all_data <- list()
all_p_scores <- list()
counter <- 1
apply_intensity <- TRUE
p_only <- FALSE
calc_mean_intensity <- TRUE
neighbor_means <- c()


k_neighbors <- 21
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

# Save the 12 wk dataframe used to calculate correlation coefs
saveRDS(test_df, paste0("Outputs/knn_analysis/data/twelve_wk_df_to_calculate_correlation_rhos_revised_with_both_n_and_p_mecp2_samples_true_k_",k_neighbors,"_neighbors", Sys.time(), ".rds"))

# Total number of iterations
num_iterations <- 1000

# Initialize the progress bar
pb <- progress_bar$new(total = num_iterations, format = "[:bar] :percent :current/:total :eta")

all_iterations_results <- list()

for (iter in 1:num_iterations) {
  
  all_p_scores <- list()
  all_data <- list()
  counter <- 1
  k <- 5
  
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
  pb$tick()
}

all_cond_data <- bind_rows(all_iterations_results)




# ----Making scatter plots with the n and p neighbors in the data included for 12 wk at various values of k-neighbors----
plot_type <- "png"
time_point <- "12wk"
fit_line <- "reg.line"
stat_used <- "spearman"
plots <- list()
k_n <- k_neighbors
counter <- 1


mean_neighbor_df_simplified <- test_df
mean_neighbor_df_simplified <- mean_neighbor_df_simplified %>%
  mutate(filename = match(filename, unique(filename)))%>%
  rename(Image = filename)



for (c in unique(mean_neighbor_df_simplified$condition)) {
  current_cond <- filter(mean_neighbor_df_simplified, condition == c)
  current_cond$Image <- as.factor(current_cond$Image)
  # current_cond <- filter(current_cond, mecp2_p == "P")
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
    ylab = "Mean Neighborhood Intensity",
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
    labs(color = "Image", fill = "Image")+
    coord_cartesian(xlim = c(0, 8000), ylim = c(0,3000))
  
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
       device = plot_type, path = "Outputs/knn_analysis/plots/",
       filename = paste0("scatter_plot_3d_combo_p_and_n_cells_included_", "_", time_point, "_",k_n, ".", plot_type),
       width = 14, height = 14, units = "in", dpi = 600, bg = "white"
)


# Saving each individual plot
for (i in 1:length(plots)) {
  ggsave(
    plot = plots[[i]],
    filename = paste0("scatter_plot_3d_p_and_n_cells_included_", i, "_", time_point, "_",k_n, ".", plot_type),
    width = 10, height = 10, units = "in", dpi = 600, bg = "white", device = plot_type, path = "Outputs/knn_analysis/plots/"
  )
}


# ----True K neighbors for 12 wk data for p only cells----
all_data <- list()
all_p_scores <- list()
counter <- 1
apply_intensity <- TRUE
p_only <- TRUE
calc_mean_intensity <- TRUE
neighbor_means <- c()


k_neighbors <- 21
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

conds <- c("NW", "NH", "SW", "SH")

for (c in conds) {
  message(paste0("Calculating rho for cond: ", c))
  current_cond <- filter(test_df, condition == c)
  cor_res <- cor.test(x = current_cond$mean, y = current_cond$positive_neighborhood_mean, method = "spearman")
  print(cor_res$estimate)
}

# Calculating rho for specific value of k. K= 3: 0.2507986 k = 5: 0.2751142 k = 7: 0.2658029, k = 9: 0.2684814, k = 11: 0.2672301, k = 21: 0.2755067 (all conds together)  
# k = 3: NW: 0.2382383 , NH: 0.1046221, SW: 0.2966918, SH: 0.2674076 / k = 5: NW: 0.2324398, NH: 0.1745932, SW: 0.3379745 , SH: 0.2584709  / k = 7: NW: 0.2042945, NH: 0.1279065, SW: 0.3592236 , SH: 0.2549389 / k = 9: NW: 0.2079402 , NH: 0.1406396, SW: 0.3504227, SH: 0.2639348 / k = 11: NW: 0.1672644, NH: 0.1659521, SW: 0.3647404, SH: 0.2748151 /
# k = 13: NW: 0.1744506 , NH: 0.168107 , SW: 0.3744558 , SH: 0.2752179   / k = 15: NW: 0.152455 , NH: 0.1442277, SW: 0.381387 , SH: 0.2910621   / k = 17: NW: 0.1530701  , NH: 0.1526787 , SW: 0.3792929, SH: 0.2736123 / k = 19: NW: 0.1289066 , NH: 0.151156 , SW: 0.3897186 , SH: 0.2787117  / k = 21: NW: 0.1304232 , NH: 0.1370063 , SW: 0.3951716, SH: 0.291192 
cor.test(x = test_df$mean, y = test_df$positive_neighborhood_mean, method = "spearman")

specific_conds_df <- data.frame(k = rep(c(3, 5, 7, 9, 11, 13, 15, 17, 19, 21), times = 4), cond = rep(c("NW", "NH", "SW", "SH"), each = 10), rho = c(0.2382383,0.2324398,
                                                                                                                                    0.2042945,0.2079402,
                                                                                                                                    0.1672644,0.1304232,
                                                                                                                                    0.1046221,0.1745932,
                                                                                                                                    0.1279065,0.1406396,
                                                                                                                                    0.1659521,0.1370063,
                                                                                                                                    0.2966918,0.3379745,
                                                                                                                                    0.3592236,0.3504227,
                                                                                                                                    0.3647404,0.3951716,
                                                                                                                                    0.2674076,0.2584709,
                                                                                                                                    0.1744506,0.152455,
                                                                                                                                    0.1530701,0.1289066,
                                                                                                                                    0.168107,0.1442277,
                                                                                                                                    0.1526787,0.151156,
                                                                                                                                    0.3744558,0.381387,
                                                                                                                                    0.3792929,0.3897186,
                                                                                                                                    0.2752179,0.2910621,
                                                                                                                                    0.2736123,0.2787117,
                                                                                                                                    0.2549389,0.2639348,
                                                                                                                                    0.2748151,0.291192))
specific_conds_df$cond <- factor(specific_conds_df$cond, levels = c("NW", "NH", "SW", "SH"))
specific_conds_df$rho <- round(specific_conds_df$rho, digits = 2)

rho_plot <- ggplot(data = specific_conds_df, aes(x = k, y = rho))+
  geom_col(fill = "steelblue")+
  facet_wrap(nrow = 2, ~cond)+
  xlab("K")+
  ylab(expression(rho))+
  ggtitle("12 Wk MeCP2+ Only")+
  theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        panel.background = element_blank(),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16, face = "bold"),
        strip.text = element_text(size = 16, face = "bold"))+
  geom_text(aes(label = rho), vjust = -0.5, size = 4)+
  coord_cartesian(ylim = c(0,0.45))


ggsave(rho_plot, device = "png", path = "Figures/", filename = "mecp2_pos_across_k.png", width = 8, height = 8, units = "in", dpi = 600, bg = "white")

# Revised 6 wk data across rho
all_data <- list()
all_p_scores <- list()
counter <- 1
apply_intensity <- TRUE
p_only <- TRUE
calc_mean_intensity <- TRUE
neighbor_means <- c()


k_neighbors <- 21
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

conds <- c("NW", "NH", "SW", "SH")

for (c in conds) {
  message(paste0("Calculating rho for cond: ", c))
  current_cond <- filter(test_df, condition == c)
  cor_res <- cor.test(x = current_cond$mean, y = current_cond$positive_neighborhood_mean, method = "spearman")
  print(cor_res$estimate)
}

# k = 3: NW: 0.3122932  , NH: 0.4806197 , SW: 0.5540698, SH: 0.5947406 / k = 5: NW: 0.2715217 , NH: 0.4437549 , SW: 0.5383963, SH: 0.5788438  / k = 7: NW: 0.2377839, NH: 0.4209924, SW: 0.528496  , SH: 0.5688027 / k = 9: NW: 0.205017  , NH: 0.4092305 , SW: 0.5152309, SH: 0.573217 / k = 11: NW: 0.1998075, NH: 0.4057064 , SW: 0.5048659 , SH: 0.5732391 /
# k = 13: NW: 0.1925582 , NH: 0.4114836 , SW: 0.501017 , SH: 0.5636041   / k = 15: NW: 0.1909383, NH: 0.4155108, SW: 0.4913128, SH: 0.5614134   / k = 17: NW: 0.1824482, NH: 0.4310985, SW: 0.4931762, SH: 0.5633391 / k = 19: NW: 0.1709586, NH: 0.4380321, SW: 0.4982802, SH: 0.5589583  / k = 21: NW: 0.1681723, NH: 0.424387 , SW: 0.4945594, SH: 0.5579393  
specific_conds_df <- data.frame(k = rep(c(3, 5, 7, 9, 11, 13, 15, 17, 19, 21), times = 4), cond = rep(c("NW", "NH", "SW", "SH"), each = 10), rho = c(0.3122932,0.2715217,
                                                                                                                                                     0.2377839,0.205017,
                                                                                                                                                     0.1998075,0.1925582,
                                                                                                                                                     0.1909383,0.1824482,
                                                                                                                                                     0.1709586,0.1681723,
                                                                                                                                                     0.4806197,0.4437549,
                                                                                                                                                     0.4209924,0.4092305 ,
                                                                                                                                                     0.4057064,0.4114836,
                                                                                                                                                     0.4155108,0.4310985,
                                                                                                                                                     0.4380321,0.424387,
                                                                                                                                                     0.5540698,0.5383963,
                                                                                                                                                     0.528496,0.5152309,
                                                                                                                                                     0.5048659,0.501017,
                                                                                                                                                     0.4913128,0.4931762,
                                                                                                                                                     0.4982802,0.4945594,
                                                                                                                                                     0.5947406,0.5788438,
                                                                                                                                                     0.5688027,0.573217,
                                                                                                                                                     0.5732391,0.5636041,
                                                                                                                                                     0.5614134,0.5633391,
                                                                                                                                                     0.5589583,0.5579393))
specific_conds_df$cond <- factor(specific_conds_df$cond, levels = c("NW", "NH", "SW", "SH"))
specific_conds_df$rho <- round(specific_conds_df$rho, digits = 2)

rho_plot <- ggplot(data = specific_conds_df, aes(x = k, y = rho))+
  geom_col(fill = "steelblue")+
  facet_wrap(nrow = 2, ~cond)+
  xlab("K")+
  ylab(expression(rho))+
  ggtitle("6 Wk MeCP2+ Only")+
  theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        panel.background = element_blank(),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16, face = "bold"),
        strip.text = element_text(size = 16, face = "bold"))+
  geom_text(aes(label = rho), vjust = -0.5, size = 4)+
  coord_cartesian(ylim = c(0,0.65))


ggsave(rho_plot, device = "png", path = "Figures/", filename = "mecp2_pos_across_k_6wk.png", width = 8, height = 8, units = "in", dpi = 600, bg = "white")









# k = 3: NW: X , NH: X, SW: X, SH: X / k = 5: NW: X, NH: X, SW: X , SH: X  / k = 7: NW: X, NH: X, SW: X , SH: X / k = 9: NW: X , NH: X, SW: X, SH: X / k = 11: NW: X, NH: X, SW: X, SH: X /
# k = 13: NW: X , NH: X , SW: X , SH: X   / k = 15: NW: X , NH: X, SW: X , SH: X   / k = 17: NW: X, NH: X , SW: X, SH: X / k = 19: NW: X , NH: X , SW: X , SH: X  / k = 21: NW: X , NH: X , SW: X, SH: X 
cor.test(x = test_df$mean, y = test_df$positive_neighborhood_mean, method = "spearman")

specific_conds_df <- data.frame(k = rep(c(3, 5, 7, 9, 11, 13, 15, 17, 19, 21), times = 4), cond = rep(c("NW", "NH", "SW", "SH"), each = 10), rho = c())
specific_conds_df$cond <- factor(specific_conds_df$cond, levels = c("NW", "NH", "SW", "SH"))
specific_conds_df$rho <- round(specific_conds_df$rho, digits = 2)

rho_plot <- ggplot(data = specific_conds_df, aes(x = k, y = rho))+
  geom_col(fill = "steelblue")+
  facet_wrap(nrow = 2, ~cond)+
  xlab("K")+
  ylab(expression(rho))+
  ggtitle("6 Wk MeCP2+ Only")+
  theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        panel.background = element_blank(),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16, face = "bold"),
        strip.text = element_text(size = 16, face = "bold"))+
  geom_text(aes(label = rho), vjust = -0.5, size = 4)+
  coord_cartesian(ylim = c(0,0.45))


ggsave(rho_plot, device = "png", path = "Figures/", filename = "mecp2_pos_across_k_6wk.png", width = 8, height = 8, units = "in", dpi = 600, bg = "white")


# Save the 12 wk dataframe used to calculate correlation coefs
saveRDS(test_df, paste0("Outputs/knn_analysis/data/twelve_wk_df_to_calculate_correlation_rhos_revised_with_p_mecp2_samples_true_k_",k_neighbors,"_neighbors", Sys.time(), ".rds"))

# Total number of iterations
num_iterations <- 1000

# Initialize the progress bar
pb <- progress_bar$new(total = num_iterations, format = "[:bar] :percent :current/:total :eta")

all_iterations_results <- list()

for (iter in 1:num_iterations) {
  
  all_p_scores <- list()
  all_data <- list()
  counter <- 1
  k <- 5
  
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
  pb$tick()
}

all_cond_data <- bind_rows(all_iterations_results)




# ----Making scatter plots with only p neighbors in the data included for 12 wk at various values of k-neighbors----
plot_type <- "png"
time_point <- "12wk"
fit_line <- "reg.line"
stat_used <- "spearman"
plots <- list()
k_n <- k_neighbors
counter <- 1


mean_neighbor_df_simplified <- test_df
mean_neighbor_df_simplified <- mean_neighbor_df_simplified %>%
  mutate(filename = match(filename, unique(filename)))%>%
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
    labs(color = "Image", fill = "Image")+
    coord_cartesian(xlim = c(0, 8000), ylim = c(0,3000))
  
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
       device = plot_type, path = "Outputs/knn_analysis/plots/",
       filename = paste0("scatter_plot_3d_combo_p_cells_only_included_", "_", time_point, "_",k_n, ".", plot_type),
       width = 14, height = 14, units = "in", dpi = 600, bg = "white"
)


# Saving each individual plot
for (i in 1:length(plots)) {
  ggsave(
    plot = plots[[i]],
    filename = paste0("scatter_plot_3d_p_cells_included_", i, "_", time_point, "_",k_n, ".", plot_type),
    width = 10, height = 10, units = "in", dpi = 600, bg = "white", device = plot_type, path = "Outputs/knn_analysis/plots/"
  )
}


# Making a barplot as we go across k for p only cells
k_neighbors <- seq(3, 30,2)
all_cors <- numeric()
apply_intensity <- TRUE
p_only <- TRUE
calc_mean_intensity <- TRUE

for (k in k_neighbors){
  for (c in conds){
    for (i in 1:length(all_imgs)) {
      all_data <- list()
      all_p_scores <- list()
      counter <- 1
      neighbor_means <- c()
      
      
      df <- all_imgs[[i]]
      df$knn_5_label <- rep(0, nrow(df))
      
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
    
    # Save the 12 wk dataframe used to calculate correlation rho
    saveRDS(test_df, paste0("Outputs/knn_analysis/data/twelve_wk_df_to_calculate_correlation_rhos_revised_with_p_mecp2_samples_true_k_",k,"_neighbors_for_barplot_across_k_", Sys.time(), ".rds"))
  }
}
