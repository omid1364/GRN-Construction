# Install and load necessary packages
install.packages("glasso")
install.packages("huge")
install.packages("glmnet")
library(glasso)
library(huge)
library(glmnet)

# Assume 'gene_expression_data' is your gene expression data frame
gene_expression_data <- as.matrix(your_data_frame)  # Convert to matrix if it's a data frame

# Standardize the data
gene_expression_data <- scale(gene_expression_data)

# Calculate the covariance matrix
cov_matrix <- cov(gene_expression_data)

# Set penalties for Lasso, Ridge, and Elastic Net
lambda <- 0.1  # Penalty parameter for Lasso and Ridge
alpha <- 0.5  # Elastic Net mixing parameter (0 = Ridge, 1 = Lasso, 0.5 = Elastic Net)

# Function to perform Glasso for different penalties
perform_glasso <- function(cov_matrix, penalty, alpha = NULL) {
  if (penalty == "lasso") {
    glasso_result <- glasso(cov_matrix, rho = lambda)
  } else if (penalty == "ridge") {
    glasso_result <- glasso(cov_matrix, rho = 0, penalize.diagonal = FALSE, zero = diag(ncol(cov_matrix)))
  } else if (penalty == "elastic_net") {
    elastic_net_result <- glmnet(cov_matrix, family = "gaussian", alpha = alpha, lambda = lambda)
    glasso_result <- glasso(cov_matrix, rho = coef(elastic_net_result)[-1])
  } else {
    stop("Invalid penalty type")
  }
  return(glasso_result$wi)  # Return the precision matrix (inverse covariance matrix)
}

# Apply the function to each penalty type
precision_matrix_lasso <- perform_glasso(cov_matrix, "lasso")
precision_matrix_ridge <- perform_glasso(cov_matrix, "ridge")
precision_matrix_elastic_net <- perform_glasso(cov_matrix, "elastic_net", alpha = alpha)

# Function to construct the network based on 95th percentile threshold
construct_network <- function(precision_matrix) {
  threshold <- quantile(abs(precision_matrix), 0.95)
  network <- abs(precision_matrix) >= threshold
  diag(network) <- 0  # Remove self-loops
  return(network)
}

# Construct networks
network_lasso <- construct_network(precision_matrix_lasso)
network_ridge <- construct_network(precision_matrix_ridge)
network_elastic_net <- construct_network(precision_matrix_elastic_net)

# Display the networks
print(network_lasso)
print(network_ridge)
print(network_elastic_net)
