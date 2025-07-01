# Step 1: Load, clean, and read the dataset for analysis (Wisconsin Breast Cancer Data)
url <- "https://archive.ics.uci.edu/ml/machine-learning-databases/breast-cancer-wisconsin/breast-cancer-wisconsin.data"
col_names <- c("ID", "Clump_Thickness", "Uniformity_of_Cell_Size", "Uniformity_of_Cell_Shape",
               "Marginal_Adhesion", "Single_Epithelial_Cell_Size", "Bare_Nuclei",
               "Bland_Chromatin", "Normal_Nucleoli", "Mitoses", "Class")
data <- read.csv(url, header = FALSE, col.names = col_names)
data$Bare_Nuclei[data$Bare_Nuclei == '?'] <- NA
data[, 2:11] <- lapply(data[, 2:11], function(x) as.numeric(as.character(x)))
data <- na.omit(data)
data <- data[, -1]

# Step 2: Load necessary libraries (install if not already installed)
if (!requireNamespace("MASS", quietly = TRUE)) install.packages("MASS")
if (!requireNamespace("mlbench", quietly = TRUE)) install.packages("mlbench")
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
if (!requireNamespace("RColorBrewer", quietly = TRUE)) install.packages("RColorBrewer")
if (!requireNamespace("gridExtra", quietly = TRUE)) install.packages("gridExtra")

library(MASS)
library(ggplot2)
library(mlbench)
library(RColorBrewer)
library(gridExtra)

# Step 3: Custom theme for plots and clustering function with plotting
my_theme <- theme_minimal() +
  theme(
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 10, hjust = 0.5),
    axis.title = element_text(size = 10, face = "bold"),
    axis.text = element_text(size = 9),
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 9),
    legend.position = "bottom",
    legend.box = "horizontal",
    panel.grid.major = element_line(color = "gray90"),
    panel.grid.minor = element_line(color = "gray95"),
    plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")
  )

fuzzy_clustering <- function(feature_name) {
  feature_data <- data[, feature_name, drop = FALSE]
  
  # Step 4: Initialize parameters and membership matrix U
  c <- 2  # Number of clusters
  m <- c(2, 3)  # Fuzziness parameters
  epsilon <- 1e-4  # Stopping criterion
  set.seed(123)
  U <- matrix(runif(nrow(feature_data) * c), nrow = nrow(feature_data), ncol = c)
  U <- U / rowSums(U)
  
  # Step 5: Functions to update cluster centers
  update_left_cluster_centers <- function(X, U, m1) {
    c <- ncol(U)
    v_L <- numeric(c)
    for (j in 1:c) {
      numerator <- sum((U[, j]^m1) * X)
      denominator <- sum(U[, j]^m1)
      v_L[j] <- numerator / denominator
    }
    return(v_L)
  }
  update_right_cluster_centers <- function(X, U, m2) {
    c <- ncol(U)
    v_R <- numeric(c)
    for (j in 1:c) {
      numerator <- sum((U[, j]^m2) * X)
      denominator <- sum(U[, j]^m2)
      v_R[j] <- numerator / denominator
    }
    return(v_R)
  }
  
  # Step 6: Calculate Euclidean distances
  calculate_distances <- function(X, centers) {
    dist <- matrix(0, nrow(X), length(centers))
    for (i in 1:nrow(X)) {
      for (j in 1:length(centers)) {
        dist[i, j] <- sum((X[i, ] - centers[j])^2)
      }
    }
    return(dist)
  }
  
  # Step 7: Calculate upper and lower membership degree for each cluster
  calculate_upper_bar_mu <- function(d_ji, d_ki, c, m1, m2) {
    ratio_sum <- sum(d_ji / d_ki)
    condition <- 1 / ratio_sum < 1 / c
    if (condition) {
      upper_bar_mu_ij <- 1 / sum((d_ji / d_ki)^(2 / (m1 - 1)))
    } else {
      upper_bar_mu_ij <- 1 / sum((d_ji / d_ki)^(2 / (m2 - 1)))
    }
    return(upper_bar_mu_ij)
  }
  calculate_lower_bar_mu <- function(d_ji, d_ki, c, m1, m2) {
    ratio_sum <- sum(d_ji / d_ki)
    condition <- 1 / ratio_sum >= 1 / c
    if (condition) {
      lower_bar_mu_ij <- 1 / sum((d_ji / d_ki)^(2 / (m1 - 1)))
    } else {
      lower_bar_mu_ij <- 1 / sum((d_ji / d_ki)^(2 / (m2 - 1)))
    }
    return(lower_bar_mu_ij)
  }
  
  # Step 8: Compute final cluster centers
  compute_final_cluster_centers <- function(v_L, v_R) {
    v <- (v_L + v_R) / 2
    return(v)
  }
  
  # Step 9: Main clustering loop
  max_iterations <- 100
  for (iteration in 1:max_iterations) {
    v_L <- update_left_cluster_centers(feature_data, U, m[1])
    v_R <- update_right_cluster_centers(feature_data, U, m[2])
    d_L <- calculate_distances(feature_data, v_L)
    d_R <- calculate_distances(feature_data, v_R)
    U_left <- matrix(0, nrow = nrow(feature_data), ncol = c)
    U_right <- matrix(0, nrow = nrow(feature_data), ncol = c)
    for (i in 1:nrow(feature_data)) {
      for (j in 1:c) {
        U_left[i, j] <- calculate_upper_bar_mu(d_L[i, j], d_L[i,], c, m[1], m[2])
        U_right[i, j] <- calculate_lower_bar_mu(d_R[i, j], d_R[i,], c, m[1], m[2])
      }
    }
    U_final <- (U_left + U_right) / 2
    v_final <- compute_final_cluster_centers(v_L, v_R)
    if (all(abs(U_final - U) < epsilon)) {
      break
    } else {
      U <- U_final
    }
  }
  
  # Step 10: Plot clusters and membership degrees
  plot_data <- data.frame(
    Feature = rep(feature_data[[1]], times = 2),
    Membership_lower = c(U_left[, 1], U_left[, 2]),
    Membership_upper = c(U_right[, 1], U_right[, 2]),
    Cluster = factor(rep(1:c, each = nrow(feature_data)))
  )
  p <- ggplot(data = plot_data, aes(x = Feature)) +
    geom_ribbon(aes(ymin = Membership_lower, ymax = Membership_upper, fill = Cluster), alpha = 0.2) +
    geom_line(aes(y = Membership_lower, color = Cluster, linetype = "Lower"), linewidth = 0.8) +
    geom_line(aes(y = Membership_upper, color = Cluster, linetype = "Upper"), linewidth = 0.8) +
    geom_vline(xintercept = v_final, linetype = "dashed", color = "#2C3E50", alpha = 0.5) +
    scale_x_continuous(limits = c(1, 10), breaks = 1:10) +
    scale_color_manual(values = c("1" = "#E74C3C", "2" = "#2ECC71"),
                       labels = c("Cluster 1", "Cluster 2")) +
    scale_fill_manual(values = c("1" = "#E74C3C", "2" = "#2ECC71"),
                      labels = c("Cluster 1", "Cluster 2")) +
    scale_linetype_manual(values = c("Lower" = "solid", "Upper" = "dashed"),
                          name = "Boundary") +
    labs(title = paste("Fuzzy Clustering Analysis:", gsub("_", " ", feature_name)),
         subtitle = "Membership Degree Boundaries and Cluster Centers",
         x = gsub("_", " ", feature_name),
         y = "Membership Degree",
         color = "Cluster",
         fill = "Cluster") +
    my_theme
  return(p)
}

# Specify features to analyze
selected_features <- c("Clump_Thickness", "Uniformity_of_Cell_Shape",
                       "Bare_Nuclei", "Single_Epithelial_Cell_Size")

# Create and print plots for selected features
plots <- list()
for(feature in selected_features) {
  plots[[feature]] <- fuzzy_clustering(feature)
}
for(feature in selected_features) {
  print(plots[[feature]])
}
