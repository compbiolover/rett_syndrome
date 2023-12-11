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




# ----Loading plotting function----
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
  ggsave(plot = p, filename = paste0(filename, ".svg"), device = "svg", width = width, height = height, units = "px", dpi = dpi, path = "~/Desktop/")
  
  # Save as PNG
  ggsave(plot = p, filename = paste0(filename, ".png"), device = "png", width = width, height = height, units = "px", dpi = dpi, path = "~/Desktop/")
  
  # Return the plot
  return(p)
}
# ----Loading function to calculate Euclidean distance between two 3D points----
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




# ----Loading original 6 wk data----
six_wk_data <- read.csv("3D_MECP2_6WO_combined_dataset.csv")
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
six_wk_data <- read.csv("Data/3D_MECP2_6WO_combined_dataset_revised.csv")
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



# ----True K=5 neighbors 6 or 12 wk data----
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

# Save the 6/12 wk dataframe used to calculate correlation coefs
saveRDS(test_df, paste0("~/Desktop/six_wk_df_to_calculate_correlation_rhos_revised_with_both_n_and_p_mecp2_samples_true_k_5_neighbors_", Sys.time(), ".rds"))

# Total number of iterations
num_iterations <- 10

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

# ----Loading 12 wk data----
twelve_wk_data <- read.csv("12WOcombined_output_MECP2.csv")
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


# ----Scatter plots with neighbors in the data included for time point----
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
       device = plot_type, path = "~/Desktop/",
       filename = paste0("scatter_plot_3d_combo_p_cells_included_", "_", time_point, ".", plot_type),
       width = 14, height = 14, units = "in", dpi = 600, bg = "white"
)


# Saving each individual plot
for (i in 1:length(plots)) {
  ggsave(
    plot = plots[[i]],
    filename = paste0("scatter_plot_3d_p_cells_included_", i, "_", time_point, ".", plot_type),
    width = 10, height = 10, units = "in", dpi = 600, bg = "white", device = plot_type, path = "~/Desktop/"
  )
}

# ----Completely random 6/12 wk KNN analysis by just condition [like initial analysis]----
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

# Save the 6/12 wk dataframe used to calculate correlation coefs
saveRDS(test_df, paste0("~/Desktop/six_wk_df_to_calculate_correlation_rhos_revised_with_only_p_mecp2_samples_", Sys.time(), ".rds"))

# Total number of iterations
num_iterations <- 10

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
num_iterations <- 10
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
x_ints <- c(0.37, 0.43, 0.50, 0.51)

# Re-order list if needed
all_results<- all_results[c("NW", "NH", "SW", "SH")]

for (cond in 1:length(conds)) {
  current_plot <- plot_permutation_test(data = all_results[[cond]], 
                                        filename = paste0(tolower(conds[cond]),"_10_random_perm_k_5_6wk_with_replacement_corrected_p_value_revised_data_p_and_n_neighbors_included_from_initial_data_finished"), 
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

ggsave(combo_plot, device = "png", path = "~/Desktop/", filename = "combo_plot_all_conditions_6wk_data_with_replacement_corrected_p_value_revised_p_and_n_neighbors_included_from_initial_data_finished.png", units = "in", width = 8, height = 8, bg = "white")
