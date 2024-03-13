# Summary of script: his code performs a comprehensive analysis of genomic data across different conditions and subregions, focusing on comparing mean values and generating heatmaps to visualize the results. The process can be summarized as follows:

# Initialization: It starts by creating a 6x6 matrix named test_mat to represent comparisons across six conditions or regions, initializing several lists for storing comparisons, estimates, significance levels, and plots. The my_mean variable is set to "mean_2", indicating the specific mean value to be analyzed.
# Comparison Preparation: The code generates a list of all possible comparisons between the conditions or regions defined in test_mat using nested for-loops.
# Data Reading and Processing: A function read_files is defined to read CSV files matching a specified pattern from a given folder, combining them into a single dataframe. This dataframe is then processed to rename certain conditions for consistency.
# Analysis Loop:
#   The script iterates over each subregion and then over each comparison within that subregion.
# For each comparison, it filters the data to include only the relevant conditions and performs a t-test to compare the mean values (mean_2) between these conditions. The p-values are adjusted using the Holm method, and significance is added based on adjusted p-values.
# Estimates, significance levels, and t-test results are stored in respective lists.
# Heatmap Generation:
#   For both overall and specific to left and right hemispheres, the code calculates differences in estimates to fill a lower triangle matrix, which is used to generate heatmaps.
# Heatmaps display differences in mean values between conditions with significant differences highlighted. Heatmaps are generated for each subregion and hemisphere, including comprehensive legends and titles to indicate the specific mean analyzed and the subregion or hemisphere in focus.
# File Output:
#   Heatmaps are saved in both PNG and SVG formats, with filenames indicating the subregion, the mean analyzed, and, where applicable, the hemisphere.
# Additionally, composite heatmaps combining multiple subregions for both left and right hemispheres are generated and saved in both formats.
# Visualization Enhancement: The script employs visual enhancements like custom color schemes, legends, and text annotations to make the heatmaps informative and easy to interpret.
# In essence, this code is designed to perform detailed comparative analysis and visualization of genomic data, focusing on how mean values differ across conditions and regions, with a special emphasis on producing high-quality, informative heatmaps for presentation and further analysis.

# Install and load devtools package and then install heatmap package
install.packages("devtools")
library(devtools)
install_github("jokergoo/ComplexHeatmap")
install.packages(c("matrixcalc", "rstatix", "tidyverse"))
library(ComplexHeatmap)
library(matrixcalc)
library(rstatix)
library(tidyverse)


# Set working directory
setwd("~/Documents/Work/Phd_program/hong_lab/Projects/rett_syndrome")

# T-tests with ComplexHeatmap
# Analyzing combined hemispheres
test_mat <- matrix(data = 0, nrow = 6, ncol = 6, dimnames = list(c("NW", "NH", "SWP1", "SWP5", "SHP1", "SHP5"), c("NW", "NH", "SWP1", "SWP5", "SHP1", "SHP5")))

comparisons <- list()
my_estimates <- list()
my_estimates_low <- list()
my_estimates_high <- list()
my_estimates_low_mod <- list()
my_sigs <- list()
my_t_tests <- list()
plots <- list()
# Update these lines with your mean of interest
my_mean <- "mean_2"
mean_title <- "Mean 2"


for (c in colnames(test_mat)) {
  for (r in rownames(test_mat)) {
    comparisons[[paste0(r, "_", c)]] <- paste0(r, "_", c)
  }
}

# Function to read in the finished files from Mclust
read_files <- function(folder_path, pattern) {
  # Get a list of all files in the folder that match the pattern
  file_list <- list.files(path = folder_path, pattern = pattern, full.names = TRUE)

  # Read in all CSV files and combine into a single dataframe
  # Assuming file_list is a vector of file names
  data <- do.call(rbind, lapply(file_list, function(x) read.csv(x, row.names = 1)))


  return(data)
}


all_plot_data <- read_files(folder_path = "Outputs/pnn_analysis/density_data/111722/", pattern = "*.csv")
all_plot_data <- all_plot_data %>%
  mutate(condition = ifelse(condition == "SH P0-P1", "SHP1", condition)) %>%
  mutate(condition = ifelse(condition == "SH P0-P5", "SHP5", condition)) %>%
  mutate(condition = ifelse(condition == "SW P0-P1", "SWP1", condition)) %>%
  mutate(condition = ifelse(condition == "SW P0-P5", "SWP5", condition))


# This code snippet is part of a larger script designed for statistical analysis and visualization of data across different subregions and conditions. Here's a summary for your collaborator:

# Subregion Iteration: The code iterates through each unique subregion within the dataset all_plot_data.
# Condition Comparisons: Within each subregion, it compares pairs of conditions specified in the comparisons list. Each comparison pair is split into two variables, cond_x and cond_y, representing the two conditions being compared.
# Data Filtering and Analysis:
# For each pair of conditions, the data is filtered to include only those entries that match either condition.
# If the filtered dataset for a comparison pair lacks entries for both conditions (i.e., less than two unique conditions are present), it assigns NA to the estimate and marks the significance as empty, indicating an inability to perform a valid comparison.
# If entries for both conditions are present, it performs a t-test to compare the mean value mean_2 between the two conditions. The p-values from the t-test are adjusted using the Holm method, and significance levels are determined based on these adjusted p-values.

# Estimate and Significance Matrix Construction:
# Estimates and significance levels from the comparisons are compiled into two matrices: estimates_mat for the estimates and sigs_mat for the significance levels, with both matrices having dimensions matching the original condition matrix test_mat.
# These matrices are then processed to reflect differences between the conditions: lower triangle values from the estimates_mat are negated and combined with upper triangle values to create a comprehensive matrix of estimate differences.
# Data Handling for Insufficient Comparisons: For comparisons with insufficient data (i.e., not both conditions represented), the code manually constructs a placeholder entry for the t-test results, ensuring consistent data structure across all comparisons.

# Final Processing:
# The code identifies and corrects zero values within the lower estimates matrix by replacing them with corresponding non-zero values from the upper estimates, ensuring a complete matrix of estimate differences for visualization or further analysis.

for (s in unique(all_plot_data$subregion)[10]) {
  current_df <- filter(all_plot_data, subregion == s)
  for (c in comparisons) {
    conds <- strsplit(c, "_")
    # cond_y <- conds[[1]][2]
    # cond_x <- conds[[1]][1]
    cond_y <- conds[[1]][1]
    cond_x <- conds[[1]][2]
    current_df_sub <- filter(current_df, condition == cond_x | condition == cond_y)
    if (length(table(current_df_sub$condition)) < 2) {
      current_comp <- c()
      current_comp$estimate <- NA
      current_comp$p.adj.signif <- ""
      my_estimates[[paste0(cond_x, cond_y)]] <- current_comp$estimate
      my_sigs[[paste0(cond_x, cond_y)]] <- current_comp$p.adj.signif
      my_t_tests[[paste0(cond_x, "_", cond_y)]] <- c(0, 0, 0, my_mean, names(table(current_df_sub$condition)), names(table(current_df_sub$condition)), unname(table(current_df_sub$condition)), unname(table(current_df_sub$condition)), 0, 1, 0, 0, 0, "T-test", "two.sided", 1, "ns", s)
    } else {
      current_comp <- NULL
      # line to change mean that is tested
      current_comp <- t_test(data = current_df_sub, formula = mean_2 ~ condition, detailed = TRUE, conf.level = 0.95, comparisons = list(c)) %>%
        adjust_pvalue(method = "holm") %>%
        add_significance("p.adj") %>%
        mutate(subregion = s)
      
      # Now doing a check to see if the t-test grouping which always does group 1 - group2 to calculate difference
      # in means has the y group in the group 1 slot and x group in the group 2 slot
      if (attributes(current_comp)$args$data[1, "condition"] != cond_y){
        print(paste0("Estimate condition: ", attributes(current_comp)$args$data[1, "condition"], " Condition y is: ", cond_y))
        print(paste0("Estimate is originally: ", current_comp$estimate))
        current_comp$estimate <- current_comp$estimate * -1
        print(paste0("Estimate is now: ", current_comp$estimate))
      }
      
      my_t_tests[[paste0(cond_x, "_", cond_y)]] <- current_comp
      my_estimates[[paste0(cond_x, "_", cond_y)]] <- current_comp$estimate
      my_sigs[[paste0(cond_x, "_", cond_y)]] <- current_comp$p.adj.signif
    }
  }

  estimates_mat <- do.call("cbind", my_estimates)
  sigs_mat <- do.call("cbind", my_sigs)
  dim(estimates_mat) <- c(6, 6)
  dim(sigs_mat) <- c(6, 6)
  dimnames(estimates_mat) <- list(rownames(test_mat), colnames(test_mat))
  dimnames(sigs_mat) <- list(rownames(test_mat), colnames(test_mat))
  lower_estimates <- lower.triangle(estimates_mat)
  upper_estimates <- upper.triangle(estimates_mat)
  lower_estimates <- lower_estimates * -1
  lower_estimates <- apply(lower_estimates, c(1, 2), as.numeric)
  indicies_to_rep <- which(lower_estimates == 0.000, arr.ind = TRUE)
  indices_to_use_in_upper <- which(upper_estimates != 0.000, arr.ind = TRUE)
  lower_estimates[indicies_to_rep] <- upper_estimates[indices_to_use_in_upper]
  
  
  # Flipping the triangles to see if that makes the heatmap correct now that the data is correct
  flip_triangles <- function(mat) {
    # Make sure mat is a matrix
    stopifnot(is.matrix(mat))
    
    # Create a copy of the matrix to preserve the original data
    new_mat <- mat
    
    # Get the upper and lower triangles, excluding the diagonal
    upper_triangle <- mat[upper.tri(mat, diag = FALSE)]
    lower_triangle <- mat[lower.tri(mat, diag = FALSE)]
    
    # Flip the upper and lower triangles
    new_mat[lower.tri(new_mat, diag = FALSE)] <- upper_triangle
    new_mat[upper.tri(new_mat, diag = FALSE)] <- lower_triangle
    
    # Return the matrix with flipped triangles
    return(new_mat)
  }
  
  # Assuming 'your_matrix' is your original matrix, apply the function
  # your_matrix <- ... # your original matrix
  # flipped_matrix <- flip_triangles(your_matrix)
  
  lower_estimates <- flip_triangles(lower_estimates)

  # This code segment is focused on generating and saving heatmaps for visualizing the differences in mean values (delta means) across various subregions and conditions, followed by creating composite images of these heatmaps. Here's a detailed summary for your collaborator:
  #
  # Heatmap Generation:
  # For each subregion, a heatmap is created using the Heatmap function from the ComplexHeatmap package in R. The heatmap visualizes the matrix of lower estimates, which represent the differences in mean values (delta means) between conditions.
  # The heatmap includes a legend, custom fonts, and title settings to enhance readability and interpretability. Specific configurations include bold fonts for titles and labels, horizontal legend orientation, and customized color for missing values (NA).
  # The heatmap's rows and columns are not clustered, ensuring the original order is maintained. Row names are placed on the left, and column names are centered and rotated for clarity. The order of the rows is explicitly set to match the predefined condition order.
  #   A custom cell function is used to annotate cells with significant differences (p < 0.05) based on t-test results, adding an additional layer of information to the heatmap.
  #   Heatmap Visualization:
  #     The heatmap is drawn with the legend positioned at the bottom. This visualization step is crucial for reviewing the heatmap before saving it to files.
  #   File Saving:
  #     The heatmap is saved in both PNG and SVG formats to ensure high-quality visual outputs suitable for different uses. The filenames reflect the subregion, analysis type (my_mean), and the format indicating a combined hemisphere view.
  #   The resolution and dimensions are set to ensure the heatmap is clearly visible and professionally presented in both digital and print formats.
  #   Composite Image Creation:
  #     After generating individual heatmaps for each subregion, the code combines these into a single composite image. This is done twice, once for saving as an SVG file and once for a PNG file, facilitating easy dissemination and presentation of the collective analysis results.
  #   The composite images are arranged in a grid layout, with the number of rows and columns specified to best display all included heatmaps. Labels and label sizes are automatically managed to maintain clarity.
  #   Purpose and Utility:
  #     The primary purpose of this code is to create detailed, visually appealing heatmaps that highlight differences in mean values across various conditions within subregions, with significant differences clearly annotated.
  #   This approach allows for an intuitive comparison of data across multiple dimensions, making it an invaluable tool for researchers looking to identify patterns, trends, and significant differences within their data.

  p1 <- Heatmap(
    matrix = lower_estimates, show_heatmap_legend = TRUE,
    heatmap_legend_param = list(
      title = "Delta Means (y-x)",
      title_gp = gpar(
        fontsize = 16,
        fontface = "bold"
      ),
      direction = "horizontal",
      title_position = "topcenter",
      labels_gp = gpar(fontsize = 14, fontface = "bold")
    ),
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    row_names_side = "left",
    column_names_centered = TRUE,
    column_names_rot = 45,
    row_order = c("SHP5", "SHP1", "SWP5", "SWP1", "NH", "NW"),
    row_title_gp = gpar(fontsize = 20, fontface = "bold"),
    row_names_gp = gpar(fontsize = 16, fontface = "bold"),
    na_col = "white", column_names_gp = gpar(
      fontsize = 16,
      fontface = "bold"
    ),
    column_title = paste0(s, " | ", mean_title),
    column_title_gp = gpar(fontsize = 22, fontface = "bold"),
    cell_fun = function(j, i, x, y, width, height, fill) {
      if (my_t_tests[[paste0(colnames(test_mat)[j], "_", rownames(test_mat)[i])]][17] < 0.05) {
        grid.text(sprintf("%s", my_t_tests[[paste0(colnames(test_mat)[j], "_", rownames(test_mat)[i])]][17]), x, y, gp = gpar(fontsize = 16))
      }
    }, rect_gp = gpar(color = "black")
  )

  p1 <- p1 %>%
    draw(heatmap_legend_side = "bottom")


  # PNG file format
  # png(file = paste0("Outputs/pnn_analysis/heatmaps/subregion_", s, "_complex_heatmap_", my_mean, "_combined_hemisphere.png"), units = "in", width = 6, height = 0.3937 * 6 + 0.7423, res = 600)
  png(file = paste0("~/Desktop/test_heatmaps/subregion_", s, "_complex_heatmap_", my_mean, "_combined_hemisphere.png"), units = "in", width = 6, height = 0.3937 * 6 + 0.7423, res = 600)
  p1 <- draw(p1, height = unit(1, "cm") * 6)
  draw(p1, height = unit(1, "cm") * 6)
  dev.off()


  # SVG file format
  # svg(file = paste0("Outputs/pnn_analysis/heatmaps/subregion_", s, "_complex_heatmap_", my_mean, "_combined_hemisphere.svg"), width = 6, height = 0.3937 * 6 + 0.7423)
  svg(file = paste0("~/Desktop/test_heatmaps/subregion_", s, "_complex_heatmap_", my_mean, "_combined_hemisphere.svg"), width = 6, height = 0.3937 * 6 + 0.7423)
  p1 <- draw(p1, height = unit(1, "cm") * 6)
  draw(p1, height = unit(1, "cm") * 6)
  dev.off()

  p1 <- grid.grabExpr(draw(p1))
  plots[[s]] <- p1
}

# SVG combo plot
combo_plot <- plot_grid(plotlist = plots, align = "hv", nrow = 5, ncol = 5, labels = "AUTO", label_size = 18)
# save_plot(filename = paste0("Outputs/pnn_analysis/heatmaps/combo_plot_complexheatmap_", my_mean, "_combined_hemisphere.svg"), plot = combo_plot, ncol = 5, nrow = 5, base_asp = 1.0, bg = "white")
save_plot(filename = paste0("~/Desktop/test_heatmaps/combo_plot_complexheatmap_", my_mean, "_combined_hemisphere.svg"), plot = combo_plot, ncol = 5, nrow = 5, base_asp = 1.0, bg = "white")

# PNG combo plot
# save_plot(filename = paste0("Outputs/pnn_analysis/heatmaps/combo_plot_complexheatmap_", my_mean, "_combined_hemisphere.png"), plot = combo_plot, ncol = 5, nrow = 5, base_asp = 1.0, bg = "white")
save_plot(filename = paste0("~/Desktop/test_heatmaps/combo_plot_complexheatmap_", my_mean, "_combined_hemisphere.png"), plot = combo_plot, ncol = 5, nrow = 5, base_asp = 1.0, bg = "white")


# This code snippet is designed to analyze and visualize data differences across various conditions within subregions and across hemispheres, focusing on statistical comparisons and heatmap generation. Here's a breakdown for your collaborator:
#
# Initialization and Setup:
# Initializes a 6x6 matrix named test_mat with zeroes, intended to represent initial data or structure for comparisons across six specified conditions or regions.
# Sets up lists for storing comparisons, estimates (along with low and high variants), significance levels, t-test results, and plots. It also defines a specific mean value (my_mean) and a title for the mean (mean_title) to be analyzed.
# Generating Comparisons:
# Populates the comparisons list with all possible pairwise combinations of conditions or regions, derived from the column and row names of test_mat.
# Data Preparation:
# Implements a function read_files to read and combine CSV files from a specified folder, matching a given pattern into a single dataframe. This dataframe is further processed to rename conditions for consistency.
# Analysis by Subregion and Hemisphere:
# Iterates through each unique subregion and hemisphere within the dataset, filtering data accordingly.
# For each pair of conditions identified earlier, it filters relevant data and assesses whether both conditions have data available. If not, placeholders are created with NA values for estimates and an indication of non-significance.
# When data for both conditions are present, performs a t-test to compare the specified mean (mean_2) between the conditions, adjusts p-values using the Holm method, and annotates significance.
# Matrix Construction for Estimates and Significance:
# Compiles results into matrices for estimates and significance, structuring these to mirror the original test_mat layout. The matrices are adjusted to only reflect differences in estimates (delta means) by modifying the lower triangle based on the comparison results.
# Fixes zero values in the lower estimates matrix by replacing them with corresponding upper estimates where applicable, ensuring a complete matrix of differences.
# Visualization Preparation:
# Prepares matrices of estimates differences for visualization. This involves negating lower triangle values and merging them with upper triangle values to accurately reflect differences across conditions.

# Analyzing the subregions by hemispheres
test_mat <- matrix(data = 0, nrow = 6, ncol = 6, dimnames = list(c("NW", "NH", "SWP1", "SWP5", "SHP1", "SHP5"), c("NW", "NH", "SWP1", "SWP5", "SHP1", "SHP5")))

comparisons <- list()
my_estimates <- list()
my_estimates_low <- list()
my_estimates_high <- list()
my_estimates_low_mod <- list()
my_sigs <- list()
my_t_tests <- list()
plots <- list()
my_mean <- "mean_2"
mean_title <- "Mean 2"


for (c in colnames(test_mat)) {
  for (r in rownames(test_mat)) {
    comparisons[[paste0(r, "_", c)]] <- paste0(r, "_", c)
  }
}

# Function to read in the finished files from Mclust
read_files <- function(folder_path, pattern) {
  # Get a list of all files in the folder that match the pattern
  file_list <- list.files(path = folder_path, pattern = pattern, full.names = TRUE)

  # Read in all CSV files and combine into a single dataframe
  data <- do.call(rbind, lapply(file_list, read.csv))

  return(data)
}


all_plot_data <- read_files(folder_path = "Outputs/pnn_analysis/density_data/111722/", pattern = "*.csv")
all_plot_data <- all_plot_data %>%
  mutate(condition = ifelse(condition == "SH P0-P1", "SHP1", condition)) %>%
  mutate(condition = ifelse(condition == "SH P0-P5", "SHP5", condition)) %>%
  mutate(condition = ifelse(condition == "SW P0-P1", "SWP1", condition)) %>%
  mutate(condition = ifelse(condition == "SW P0-P5", "SWP5", condition))

lh_plots <- list()
rh_plots <- list()


for (s in unique(all_plot_data$subregion)) {
  current_subregion <- filter(all_plot_data, subregion == s)
  for (h in unique(all_plot_data$hemisphere)) {
    current_df <- filter(current_subregion, hemisphere == h)
    for (c in comparisons) {
      conds <- strsplit(c, "_")
      cond_y <- conds[[1]][2]
      cond_x <- conds[[1]][1]
      current_df_sub <- filter(current_df, condition == cond_x | condition == cond_y)
      print(head(current_df_sub))
      if (length(table(current_df_sub$condition)) < 2) {
        current_comp <- c()
        current_comp$estimate <- NA
        current_comp$p.adj.signif <- ""
        my_estimates[[paste0(cond_x, cond_y)]] <- current_comp$estimate
        my_sigs[[paste0(cond_x, cond_y)]] <- current_comp$p.adj.signif
        my_t_tests[[paste0(cond_x, "_", cond_y)]] <- c(0, 0, 0, my_mean, names(table(current_df_sub$condition)), names(table(current_df_sub$condition)), unname(table(current_df_sub$condition)), unname(table(current_df_sub$condition)), 0, 1, 0, 0, 0, "T-test", "two.sided", 1, "ns")
      } else {
        current_comp <- NULL
        # line to change mean that is tested
        current_comp <- t_test(data = current_df_sub, formula = mean_2 ~ condition, detailed = TRUE, conf.level = 0.95) %>%
          adjust_pvalue(method = "holm") %>%
          add_significance("p.adj")
        my_t_tests[[paste0(cond_x, "_", cond_y)]] <- current_comp
        my_estimates[[paste0(cond_x, "_", cond_y)]] <- current_comp$estimate
        my_sigs[[paste0(cond_x, "_", cond_y)]] <- current_comp$p.adj.signif
      }
    }

    estimates_mat <- do.call("cbind", my_estimates)
    sigs_mat <- do.call("cbind", my_sigs)
    dim(estimates_mat) <- c(6, 6)
    dim(sigs_mat) <- c(6, 6)
    dimnames(estimates_mat) <- list(rownames(test_mat), colnames(test_mat))
    dimnames(sigs_mat) <- list(rownames(test_mat), colnames(test_mat))
    lower_estimates <- lower.triangle(estimates_mat)
    upper_estimates <- upper.triangle(estimates_mat)
    lower_estimates <- lower_estimates * -1
    lower_estimates <- apply(lower_estimates, c(1, 2), as.numeric)
    indicies_to_rep <- which(lower_estimates == 0.000, arr.ind = TRUE)
    indices_to_use_in_upper <- which(upper_estimates != 0.000, arr.ind = TRUE)
    lower_estimates[indicies_to_rep] <- upper_estimates[indices_to_use_in_upper]


    # Heatmap like before

    p1 <- Heatmap(
      matrix = lower_estimates, show_heatmap_legend = TRUE,
      heatmap_legend_param = list(
        title = "Delta Means (y-x)",
        title_gp = gpar(
          fontsize = 16,
          fontface = "bold"
        ),
        direction = "horizontal",
        title_position = "topcenter",
        labels_gp = gpar(fontsize = 14, fontface = "bold")
      ),
      cluster_rows = FALSE,
      cluster_columns = FALSE,
      row_names_side = "left",
      column_names_centered = TRUE,
      column_names_rot = 45,
      row_order = c("SHP5", "SHP1", "SWP5", "SWP1", "NH", "NW"),
      row_title_gp = gpar(fontsize = 20, fontface = "bold"),
      row_names_gp = gpar(fontsize = 16, fontface = "bold"),
      na_col = "white", column_names_gp = gpar(
        fontsize = 16,
        fontface = "bold"
      ),
      column_title = paste0(s, " | ", mean_title, " | ", h),
      column_title_gp = gpar(fontsize = 22, fontface = "bold"),
      cell_fun = function(j, i, x, y, width, height, fill) {
        if (my_t_tests[[paste0(colnames(test_mat)[j], "_", rownames(test_mat)[i])]][17] < 0.05) {
          grid.text(sprintf("%s", my_t_tests[[paste0(colnames(test_mat)[j], "_", rownames(test_mat)[i])]][17]), x, y, gp = gpar(fontsize = 16))
        }
      }, rect_gp = gpar(color = "black")
    )

    p1 <- p1 %>%
      draw(heatmap_legend_side = "bottom")

    # PNG file format
    png(file = paste0("Outputs/pnn_analysis/heatmaps/subregion_", s, "_complex_heatmap_", my_mean, "_", h, ".png"), units = "in", width = 6, height = 0.3937 * 6 + 0.7423, res = 600)
    p1 <- draw(p1, height = unit(1, "cm") * 6)
    draw(p1, height = unit(1, "cm") * 6)
    dev.off()


    # SVG file format
    svg(file = paste0("Outputs/pnn_analysis/heatmaps/subregion_", s, "_complex_heatmap_", my_mean, "_", h, ".svg"), width = 6, height = 0.3937 * 6 + 0.7423)
    p1 <- draw(p1, height = unit(1, "cm") * 6)
    draw(p1, height = unit(1, "cm") * 6)
    dev.off()



    p1 <- grid.grabExpr(draw(p1))
    if (h == "LH") {
      lh_plots[[s]] <- p1
    } else {
      rh_plots[[s]] <- p1
    }
  }
}

for (h in unique(all_plot_data$hemisphere)) {
  if (h == "LH") {
    combo_plot <- plot_grid(plotlist = lh_plots, align = "hv", nrow = 5, ncol = 5, labels = "AUTO", label_size = 18)
    save_plot(filename = paste0("Outputs/pnn_analysis/heatmaps/combo_plot_complexheatmap_", my_mean, "_", h, ".svg"), plot = combo_plot, ncol = 5, nrow = 5, base_asp = 1.0, bg = "white")
    save_plot(filename = paste0("Outputs/pnn_analysis/heatmaps/combo_plot_complexheatmap_", my_mean, "_", h, ".png"), plot = combo_plot, ncol = 5, nrow = 5, base_asp = 1.0, bg = "white")
  } else {
    combo_plot <- plot_grid(plotlist = rh_plots, align = "hv", nrow = 5, ncol = 5, labels = "AUTO", label_size = 18)
    save_plot(filename = paste0("Outputs/pnn_analysis/heatmaps/combo_plot_complexheatmap_", my_mean, "_", h, ".svg"), plot = combo_plot, ncol = 5, nrow = 5, base_asp = 1.0, bg = "white")
    save_plot(filename = paste0("Outputs/pnn_analysis/heatmaps/combo_plot_complexheatmap_", my_mean, "_", h, ".png"), plot = combo_plot, ncol = 5, nrow = 5, base_asp = 1.0, bg = "white")
  }
}
