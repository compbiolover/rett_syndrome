library(ComplexHeatmap)
library(matrixcalc)
library(rstatix)
library(tidyverse)

test_mat <- matrix(data = 0, nrow = 6, ncol = 6, dimnames = list(c("NW", "NH", "SWD2", "SWD6", "SHD2", "SHD6"), c("NW", "NH", "SWD2", "SWD6", "SHD2", "SHD6")))

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

# Generate possible combinations for pairwise comparisons later
for (c in colnames(test_mat)) {
  for (r in rownames(test_mat)) {
    comparisons[[paste0(r, "_", c)]] <- paste0(r, "_", c)
  }
}

# Jacob's new data
all_plot_data <- read.csv("Data/combined_JE_results.csv", row.names = 1)
all_plot_data <- all_plot_data %>% select(-X)

# Remove cohort that has only one condition for all samples (SHD2 condition)
all_plot_data <- filter(all_plot_data, cohort != "072123X")


for (s in unique(all_plot_data$subregion)[7]) {
  current_df <- filter(all_plot_data, subregion == s)
  for (co in unique(all_plot_data$cohort)[2]) {
    current_df_coh <- filter(current_df, cohort == co)
    print(paste0("On cohort ", co))

    for (c in comparisons) {
      conds <- strsplit(c, "_")[[1]]
      cond_y <- conds[1]
      cond_x <- conds[2]

      current_df_sub <- filter(current_df_coh, condition == cond_x | condition == cond_y)

      if (cond_x != cond_y) {
        condition_table <- table(current_df_sub$condition)

        # Check if any of the two conditions have zero observations or they have only one condition (self comparisons or missing condition)
        if (any(condition_table == 0) || length(condition_table) < 2) {
          # Skip t-test if either condition has zero observations
          my_estimates[[paste0(cond_x, "_", cond_y)]] <- NA
          my_sigs[[paste0(cond_x, "_", cond_y)]] <- NA
          my_t_tests[[paste0(cond_x, "_", cond_y)]] <- as_tibble(list(
            estimate = NA,
            estimate1 = NA,
            estimate2 = NA,
            .y. = NA,
            group1 = NA,
            group2 = NA,
            n1 = NA,
            n2 = NA,
            statistic = NA,
            p = NA,
            df = NA,
            conf.low = NA,
            conf.high = NA,
            method = NA,
            alternative = NA,
            p.adj = NA,
            p.adj.signif = NA,
            subregion = NA
          ))
        } else {
          # Proceed with t-test
          current_comp <- t_test(data = current_df_sub, formula = mean_2 ~ condition, detailed = TRUE, conf.level = 0.95) %>%
            adjust_pvalue(method = "holm") %>%
            add_significance("p.adj") %>%
            mutate(subregion = s)

          # Checking one final time if the group order in t-test results matches the expected order (cond_y - cond_x)
          group1 <- current_comp$group1 # Get the name of the first group from t_test results
          group2 <- current_comp$group2 # Get the name of the second group from t_test results

          # If the group order is not the expected one, reverse the estimate sign and the data related to each group
          if (group1 != cond_y | group2 != cond_x) {
            print(paste0("T-test group1 condition: ", group1, ", T-test group2 condition: ", group2))
            print(paste0("Estimate is originally: ", current_comp$estimate))
            print(paste0("Our initial pairing is ", cond_y, "_", cond_x))
            group1 <- cond_y
            group2 <- cond_x
            current_comp$estimate <- current_comp$estimate * -1
            current_comp$group1 <- cond_y
            current_comp$group2 <- cond_x
            estimate2 <- current_comp$estimate1
            estimate1 <- current_comp$estimate2
            n2 <- current_comp$n1
            n1 <- current_comp$n2
            current_comp$estimate1 <- estimate1
            current_comp$estimate2 <- estimate2
            current_comp$n1 <- n1
            current_comp$n2 <- n2
            print(paste0("After flipping the sign to match our matrix the estimate is: ", current_comp$estimate))
            print(paste0("We have updated the groups and data of the my_t_tests opbject to reflect the correct grouping. Group1 is now: ", group1, " and Group2 is: ", group2))
          }

          # Store results in lists that follow earlier convention (y - x)
          my_t_tests[[paste0(cond_y, "_", cond_x)]] <- current_comp
          my_estimates[[paste0(cond_y, "_", cond_x)]] <- current_comp$estimate
          my_sigs[[paste0(cond_y, "_", cond_x)]] <- current_comp$p.adj.signif
        }
      } else {
        # Skip t-test if conditions are the same
        my_estimates[[paste0(cond_y, "_", cond_x)]] <- NA
        my_sigs[[paste0(cond_y, "_", cond_x)]] <- NA
        my_t_tests[[paste0(cond_y, "_", cond_x)]] <- as_tibble(list(
          estimate = NA,
          estimate1 = NA,
          estimate2 = NA,
          .y. = NA,
          group1 = NA,
          group2 = NA,
          n1 = NA,
          n2 = NA,
          statistic = NA,
          p = NA,
          df = NA,
          conf.low = NA,
          conf.high = NA,
          method = NA,
          alternative = NA,
          p.adj = NA,
          p.adj.signif = NA,
          subregion = NA
        ))
      }
    }

    # Now processing all of the data for heatmap
    estimates_mat <- do.call("cbind", my_estimates)
    sigs_mat <- do.call("cbind", my_sigs)
    dim(estimates_mat) <- c(6, 6)
    dim(sigs_mat) <- c(6, 6)
    dimnames(estimates_mat) <- list(rownames(test_mat), colnames(test_mat))
    dimnames(sigs_mat) <- list(rownames(test_mat), colnames(test_mat))
    upper_estimates <- upper.triangle(estimates_mat)
    
    # Function to mirror the upper triangle to the lower triangle and flip the sign so it can be read in proper direction
    mirror_upper_to_lower <- function(matrix) {
      # Ensure the matrix is square
      if (!all(dim(matrix)[1] == dim(matrix)[2])) {
        stop("The matrix must be square.")
      }
      
      # Copy the upper triangle to the lower triangle
      matrix[lower.tri(matrix)] <- -1*(t(matrix)[lower.tri(matrix)])
      
      return(matrix)
    }
    
    lower_estimates <- mirror_upper_to_lower(upper_estimates)
    

    # Heatmap for Jacob's new data (by cohort)
    p1 <- Heatmap(
      matrix = lower_estimates, show_heatmap_legend = TRUE,
      heatmap_legend_param = list(
        title = expression(bold(Delta * " Means (y - x)")),
        title_gp = gpar(
          fontsize = 16,
          fontface = "bold"
        ),
        direction = "horizontal",
        title_position = "topcenter",
        width = unit(0.8, "snpc"),
        height = unit(0.4, "snpc"),
        labels_gp = gpar(fontsize = 14, fontface = "bold")
      ),
      cluster_rows = FALSE,
      cluster_columns = FALSE,
      column_names_centered = TRUE,
      column_names_side = "bottom",
      column_names_rot = 45,
      column_order = c("NW", "NH", "SWD2", "SWD6", "SHD2", "SHD6"),
      row_order = c("SHD6", "SHD2", "SWD6", "SWD2", "NH", "NW"),
      row_title_gp = gpar(fontsize = 20, fontface = "bold"),
      row_names_gp = gpar(fontsize = 16, fontface = "bold"),
      row_names_side = "left",
      na_col = "white", column_names_gp = gpar(
        fontsize = 16,
        fontface = "bold"
      ),
      column_title = paste0(s, " | ", mean_title, " | ", co),
      column_title_gp = gpar(fontsize = 22, fontface = "bold"),
      cell_fun = function(j, i, x, y, width, height, fill) {
        key <- paste0(rownames(test_mat)[i], "_", colnames(test_mat)[j])
        if (my_t_tests[[key]][16] < 0.05 && !is.na(my_t_tests[[key]][16])) {
          sig_stars <- my_t_tests[[key]][17]
          grid.text(sig_stars, x, y, gp = gpar(fontsize = 16))
        }
      },
      rect_gp = gpar(color = "black")
    )

    p1 <- p1 %>%
      draw(heatmap_legend_side = "bottom")


    # PNG file format
    png(file = paste0("~/Desktop/test_heatmaps/subregion_", s, "_complex_heatmap_", my_mean, "_combined_hemisphere_jacob_data_cohort_", co, ".png"), units = "in", width = 6, height = 0.3937 * 6 + 0.7423, res = 600)
    p1 <- draw(p1, height = unit(1, "cm") * 6)
    draw(p1, height = unit(1, "cm") * 6)
    dev.off()


    # SVG file format
    svg(file = paste0("~/Desktop/test_heatmaps/subregion_", s, "_complex_heatmap_", my_mean, "_combined_hemisphere_jacob_data_cohort_", co, ".svg"), width = 6, height = 0.3937 * 6 + 0.7423)
    p1 <- draw(p1, height = unit(1, "cm") * 6)
    draw(p1, height = unit(1, "cm") * 6)
    dev.off()
  }
}

