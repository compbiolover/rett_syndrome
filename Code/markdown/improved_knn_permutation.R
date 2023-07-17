# Name: improved_knn_permutation.R
# Purpose: Re-write of both test_3d_knn.R and parts of knn_analysis_2.Rmd. To make more readable, reproducible, and
# modular code. To also speed up code
# Author: Andrew Willems <awillems@vols.utk.edu>

# Loading packages
pacman::p_load(
  car, # Tools for performing specific statistical tests
  caret, # Tools for training and evaluating predictive models
  class, # Functions for nearest neighbor classification
  deldir, # Delaunay triangulation and Dirichlet (Voronoi) tesselation
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


#' Euclidean Distance 2D Function
#'
#' Calculate the Euclidean distance between two points in 2 dimensions (x, y).
#'
#' @param x1, y1 Numeric vectors representing the coordinates of the first point.
#' @param x2, y2 Numeric vectors representing the coordinates of the second point.
#'
#' @return Numeric scalar representing the Euclidean distance between the two points.
#'
#' @examples
#' euclidean_distance(x1 = 10, y1 = 20, x2 = 12, y2 = 22)
#' Answer is: 2.828427
#' @references
#' For more information on Euclidean distance, see:
#' \url{https://en.wikipedia.org/wiki/Euclidean_distance}
#'
#' @keywords math distance euclidean 2d
#' @export
euclidean_distance <- function(x1, y1, x2, y2) {
  # Check if inputs are numeric vectors
  if (!is.numeric(x1) || !is.numeric(y1) || !is.numeric(x2) || !is.numeric(y2)) {
    stop("Inputs must be numeric vectors")
  }
  
  # Check if inputs have the same length
  if (length(x1) == length(y1) && length(x2) == length(y2)) {
    sqrt(sum((x2 - x1)^2 + (y2 - y1)^2))
  } else {
    stop("Vectors must be of the same length")
  }
}


#' Calculate Euclidean Distance in 3D Space
#'
#' Calculate the Euclidean distance between two points in 3D space.
#'
#' @param x1, y1, z1 Numeric vectors representing the coordinates of the first point.
#' @param x2, y2, z2 Numeric vectors representing the coordinates of the second point.
#'
#' @return Numeric vector representing the Euclidean distance between the corresponding points.
#'
#' @examples
#' euclidean_distance_3d(x1 = 10, y1 = 20, z1 = 5, x2 = 12, y2 = 22, z2 = 3)
#' Answer is: 3.464102
#' @references
#' For more information on Euclidean distance, see:
#' \url{https://en.wikipedia.org/wiki/Euclidean_distance}
#'
#' @keywords math distance euclidean 3d
#' @export
euclidean_distance_3d <- function(x1, y1, z1, x2, y2, z2) {
  # Check if inputs are numeric vectors
  if (!is.numeric(x1) || !is.numeric(y1) || !is.numeric(z1) ||
      !is.numeric(x2) || !is.numeric(y2) || !is.numeric(z2)) {
    stop("Inputs must be numeric vectors")
  }
  sqrt((x2 - x1)^2 + (y2 - y1)^2 + (z2 - z1)^2)
}



# Permutation analysis
permutation_knn_analysis <- function(all_imgs,
                                     num_nearest = 30,
                                     custom_seed = 1689014354,
                                     k = 5,
                                     apply_intensity = TRUE,
                                     p_only = TRUE,
                                     calc_mean_intensity = TRUE,
                                     num_iterations = 1000,
                                     num_cores = NULL,
                                     random_num_neighbors = FALSE,
                                     max_randomness = FALSE) {
  
  # Create an empty vector to store the correlation values
  all_rhos_stored <- vector()
  
  # Create an empty vector to store the p-values
  all_p_values_stored <- vector()
  
  # Set the random seed if provided
  if (!is.null(custom_seed)) {
    set.seed(custom_seed)
    cat("Custom random seed:", custom_seed, "\n")
  } else {
    # Generate a random seed using the system time
    random_seed <- as.integer(Sys.time())
    set.seed(random_seed)
    cat("Generated random seed:", random_seed, "\n")
  }
  
  # Determine the number of CPU cores to use
  if (is.null(num_cores)) {
    num_cores <- parallel::detectCores() - 1
    cat("Number of cores used:", num_cores, "\n")
  } else {
    num_cores <- num_cores
    cat("Number of cores used:", num_cores, "\n")
  }
  
  
  
  # Create a progress bar
  pb <- progress_bar$new(total = num_iterations, width = 60)
  
  # List to store results from the permutation
  results <- vector("list", num_iterations)
  
  # Perform the permutation analysis
  img_counter <- 1
  result_index <- 1
  for (iteration in 1:num_iterations) {
    # Create an empty vector to store the iteration results
    iteration_rhos <- vector()
    iteration_p_values <- vector()
    
    # Loop through each image
    for (i in seq_along(all_imgs)) {
      df <- all_imgs[[i]]
      n <- nrow(df)
      
      all_p_scores <- list()
      all_data <- list()
      counter <- 1
      rhos <- numeric(length(all_imgs))
      p_values <- numeric(length(all_imgs))
      
      # For each row within an image
      for (r in seq_len(n)) {
        observation <- df[r, ]
        other_points <- df[-r, ]
        
        if (random_num_neighbors) {
          # First calculate distances
          distances <- euclidean_distance_3d(
            observation$x, observation$y, observation$z,
            other_points$x, other_points$y, other_points$z
          )
          # Then subset to the neighbors within num_nearest with order
          k_nearest_indices <- order(distances)[1:num_nearest]
          # Then sample from these num_nearest neighbors randomly
          k_nearest_indices <- sample(k_nearest_indices, replace = TRUE)
          # Then pass the randomly selected indices from the num_nearest neighbors to the other points
          k_nearest_points <- other_points[k_nearest_indices, ]
          
          # For max randomness
          if (max_randomness){
            k_nearest_indices <- sample(nrow(other_points), k, replace = TRUE)
            k_nearest_points <- other_points[k_nearest_indices, ]
          }
        } else {
          # Perform regular KNN analysis with no randomness
          distances <- euclidean_distance_3d(
            observation$x, observation$y, observation$z,
            other_points$x, other_points$y, other_points$z
          )
          k_nearest_indices <- order(distances)[1:k]
          k_nearest_points <- other_points[k_nearest_indices, ]
        }
        
        # Calculating other parts of the score beyond just KNN
        add_p_score <- sum((k_nearest_points$mecp2_p == "P")) / k
        observation$knn_label <- add_p_score
        all_p_scores[[counter]] <- add_p_score
        
        # Calculate distances with intensity information
        if (apply_intensity) {
          distances <- distances[k_nearest_indices]
          distances_updated <- rescale(distances * k_nearest_points$mean, to = c(0, 1))
          distances_updated <- sort(distances_updated, decreasing = FALSE)
          
          # Calculate updated distances with intensity information for just P neighbors
          if (p_only) {
            k_nearest_points <- k_nearest_points[k_nearest_points$mecp2_p == "P", ]
            if (nrow(k_nearest_points) == 0) {
              writeLines("This sample has only N neighbors.\nIt is being excluded from this analysis.")
              next
            }
          }
          if (calc_mean_intensity) {
            mean_of_neighbors <- mean(k_nearest_points$mean)
            observation$positive_neighborhood_mean <- mean_of_neighbors
            all_data <- c(all_data, list(observation))
            counter <- counter + 1
          }
        }
      }
      
      # Binding all observations of an image together
      df <- bind_rows(all_data)
      
      # Performing the correlation test between the mean and the positive (mecp2_p = P) neighborhood mean of
      # all rows in a particular image
      cor_res <- cor.test(
        x = df$mean,
        y = df$positive_neighborhood_mean,
        method = "spearman"
      )
      rhos[i] <- as.numeric(as.vector(unlist(cor_res$estimate)))
      p_values[i] <- cor_res$p.value
      
      # Store the iteration results in the results list
      results[[paste0("Iteration ", result_index)]] <- data.frame(image = img_counter,
                                                              rho = rhos[i],
                                                              p_value = p_values[i],
                                                              condition = unique(df$condition),
                                                              time = unique(df$time),
                                                              k = k,
                                                              k_search_space = num_nearest,
                                                              num_iterations = num_iterations)
      
      # img_counter <- img_counter + 1
      result_index <- result_index + 1
      
      
    }
    
    # iteration_rhos <- c(iteration_rhos, rhos)
    # iteration_p_values <- c(iteration_p_values, p_values)
    
    # Update the progress bar
    pb$tick()
    img_counter <- img_counter + 1
    # print(results)
  }
  
  # Combine the iteration results into a single data frame
  result_df <- do.call(rbind, results)
  
  # Return the result data frame
  return(result_df)
}







# 6 week data
six_wk_data <- read.csv("Data/Old data/Jacob/3D_MECP2_6WO_combined_dataset.csv")
six_wk_data <- six_wk_data %>%
  select(Filename, Mean, CX..pix., CY..pix., CZ..pix., MECP2) %>%
  rename_all(tolower) %>%
  rename(
    x = cx..pix.,
    y = cy..pix.,
    z = cz..pix.,
    mecp2_p = mecp2
  ) %>%
  mutate(hemisphere = str_extract(filename, "(LH|RH)")) %>%
  mutate(condition = str_extract(filename, "(NW|NH|SW|SH)")) %>%
  filter(mecp2_p == "P") %>%
  mutate(id = 1:nrow(.)) %>%
  mutate(time = 6) %>%
  # Modify this line to filter by different condition
  filter(condition == "SW")

six_wk_data <- six_wk_data %>% group_by(filename)
all_imgs <- six_wk_data %>% group_split(six_wk_data)

# 12 week data
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
  # Modify this line to filter by different condition
  filter(condition == "NW")

twelve_wk_data <- twelve_wk_data %>% group_by(image)
all_imgs <- twelve_wk_data %>% group_split(twelve_wk_data)


# NW 10000 permutations at 30 search space
nw_results_10000_30 <- permutation_knn_analysis(all_imgs,
                                         num_nearest = 30,
                                         custom_seed = 1689016866,
                                         k = 5,
                                         num_iterations = 10000,
                                         num_cores = 1,
                                         random_num_neighbors = TRUE,
                                         max_randomness = FALSE,
                                         apply_intensity = TRUE,
                                         p_only = TRUE,
                                         calc_mean_intensity = TRUE)


# NW 10000 permutations at 5-50 search space with steps of 5
search_space <- seq(5, 50, 5)
all_nw_searches <- list()
for (s in search_space){
  current_search <- permutation_knn_analysis(all_imgs,
                                                  num_nearest = s,
                                                  custom_seed = 1689016866,
                                                  k = 5,
                                                  num_iterations = 10000,
                                                  num_cores = 1,
                                                  random_num_neighbors = TRUE,
                                                  max_randomness = FALSE,
                                                  apply_intensity = TRUE,
                                                  p_only = TRUE,
                                                  calc_mean_intensity = TRUE) 
  
  # Append the result of the new search space to all results
  all_nw_searches[[paste0("Search ", s)]] <- current_search
}




# NH 10000 permutations at 5-50 search space with steps of 5
search_space <- seq(5, 50, 5)
all_nh_searches <- list()
for (s in search_space){
  current_search <- permutation_knn_analysis(all_imgs,
                                             num_nearest = s,
                                             custom_seed = 1689016866,
                                             k = 5,
                                             num_iterations = 10000,
                                             num_cores = 1,
                                             random_num_neighbors = TRUE,
                                             max_randomness = FALSE,
                                             apply_intensity = TRUE,
                                             p_only = TRUE,
                                             calc_mean_intensity = TRUE) 
  
  # Append the result of the new search space to all results
  all_nh_searches[[paste0("Search ", s)]] <- current_search
}


# SW 10000 permutations at 5-50 search space with steps of 5
search_space <- seq(5, 50, 5)
all_sw_searches <- list()
for (s in search_space){
  current_search <- permutation_knn_analysis(all_imgs,
                                             num_nearest = s,
                                             custom_seed = 1689016866,
                                             k = 5,
                                             num_iterations = 1000,
                                             num_cores = 1,
                                             random_num_neighbors = TRUE,
                                             max_randomness = FALSE,
                                             apply_intensity = TRUE,
                                             p_only = TRUE,
                                             calc_mean_intensity = TRUE) 
  
  # Append the result of the new search space to all results
  all_sw_searches[[paste0("Search ", s)]] <- current_search
}





#' Plot histogram for permutation test results
#'
#' This function creates a histogram plot based on the provided data, showing the results of a permutation test.
#'
#' @param data The dataframe containing the data to be plotted.
#' @param filename The filename for saving the plot (without extension).
#' @param x_intercept The x-coordinate value for the vertical line to be added to the plot.
#' @param plot_size The desired size of the plot (in inches).
#' @param base_plot_size The base plot size to be used for scaling text elements.
#' @param dpi The DPI (dots per inch) for the plot.
#' @param hist_color The color of the histogram bars.
#' @param hist_fill The fill color of the histogram bars.
#' @param num_bins The number of bins to use in the histogram.
#' @param line_color The color of the vertical line.
#' @param text_color The color of the annotation text.
#'
#' @return The created ggplot object representing the histogram plot.
#'
#' @import ggplot2
#' @import grid
#' @import scales
#'
#' @export
plot_permutation_test <- function(data, filename, x_intercept, plot_size = 10, base_plot_size = 8, dpi = 600, hist_color = "steelblue",
                                  hist_fill = "steelblue", num_bins = 30, line_color = "red", text_color = "red") {
  # Required packages
  library(ggplot2)
  library(grid)
  library(scales)

  time_point <- unique(data$time)
  k_value <- unique(data$k)
  num_permutations <- unique(data$num_iterations)
  num_ks <- unique(data$k_search_space)
  
  
  # Calculate the scaling factor for text elements based on plot size and base plot size
  text_scale <- plot_size / base_plot_size

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
      axis.line = element_line(colour = "black", linewidth = 1.25 * text_scale)
    ) +
    ggtitle(paste0("K = ", k_value, " | ", unique(data$condition), " | ", num_permutations, " Permutations | \n", time_point, " WK | ", num_ks, " KNNs Searched" )) +
    xlab(expression(rho)) +
    ylab("Count") +
    scale_y_continuous(expand = c(0, 0)) +
    geom_vline(xintercept = x_intercept, color = line_color)


  # Calculate the width and height for the plot
  width <- plot_size * dpi
  height <- plot_size * dpi

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

# Call the function for each dataset with custom parameters
nh_plot <- plot_permutation_test(nh_data, "nh_random_knn_vs_rho_histogram", 0.17,
  plot_size = 4, base_plot_size = 8, dpi = 300, hist_color = "blue",
  hist_fill = "lightblue", line_color = "green", text_color = "green"
)
nw_plot <- plot_permutation_test(all_nw_searches[[1]], "nw_random_knn_vs_rho_histogram_10_10000", 0.27,
  plot_size = 2.5, base_plot_size = 8, dpi = 600, hist_color = "blue",
  hist_fill = "lightblue", line_color = "red", text_color = "red"
)


# Looping through all outputs again to make comboplot for easier comparison of results
combo_plot_list <- vector("list", length = length(all_nw_searches))
for (s in 1:length(all_nw_searches)){
  combo_plot_list[[s]] <- plot_permutation_test(all_nw_searches[[s]], paste0("nw_random_knn_vs_rho_histogram_",s,"_10000"), 0.27,
                                                plot_size = 2.5, base_plot_size = 8, dpi = 600, hist_color = "blue",
                                                hist_fill = "lightblue", line_color = "red", text_color = "red"
  )
}

# Making combo plot
combo_plot <- ggarrange(plotlist = combo_plot_list,
                        ncol = 2,
                        nrow = 5,
                        labels = "AUTO",
                        align = "hv", 
                        font.label = list(size = 18, face = "bold"))


# Saving combo plot
ggsave(plot = combo_plot, device = "png", path = "Outputs/knn_analysis/plots/", filename = "nw_random_searches_10000_combo.png", units = "in", dpi = 600, width = 8, height = 8)



sh_plot <- plot_permutation_test(sh_data, "sh_random_knn_vs_rho_histogram_test", 0.58,
  plot_size = 2.5, base_plot_size = 8, dpi = 300, hist_color = "gray",
  hist_fill = "lightgray", line_color = "brown", text_color = "brown", num_bins = 5
)
