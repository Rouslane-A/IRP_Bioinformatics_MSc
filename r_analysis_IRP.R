library(TCGAbiolinks)
library(ggplot2)
library(gplots)
library(Boruta)
library(superml)
library(xgboost)
library(caret)
library(randomForest)
library(glmnet)
library(FactoMineR)
library(factoextra)
library(cluster)
library(tidyverse)
library(msigdbr)
library(clusterProfiler)
library(enrichplot)
library(ggupset)
library(DOSE)
library(pROC)
library(e1071)
library(maftools)
library(RColorBrewer)
library(pheatmap)

set.seed(42)

# Data Acquisition
#projects <- TCGAbiolinks:::getGDCprojects()
#print(projects)

# Download RNA-seq data for primary tumors and solid tissue normal samples

# primary tumours
query_rnaseq_primary_tumours <- GDCquery(project = "TCGA-HNSC",
                                         data.category = 'Transcriptome Profiling',
                                         sample.type = 'Primary Tumor',
                                         data.format = 'TSV',
                                         data.type = 'Gene Expression Quantification',
                                         workflow.type = 'STAR - Counts')

GDCdownload(query_rnaseq_primary_tumours, method = "api", directory = "GDCdata_TCGA_HNSC_RNASEQ_primary_tumours")

hnsc_rnaseq_tp <- GDCprepare(query_rnaseq_primary_tumours, directory = "GDCdata_TCGA_HNSC_RNASEQ_primary_tumours")

# solid tissue normal
query_rnaseq_solid_tissue_normal <- GDCquery(project = "TCGA-HNSC",
                                             data.category = 'Transcriptome Profiling',
                                             sample.type = 'Solid Tissue Normal',
                                             data.format = 'TSV',
                                             data.type = 'Gene Expression Quantification',
                                             workflow.type = 'STAR - Counts')

GDCdownload(query_rnaseq_solid_tissue_normal, method = "api", directory = "GDCdata_TCGA_HNSC_RNASEQ_solid_tissue_normal")

hnsc_rnaseq_nt <- GDCprepare(query_rnaseq_solid_tissue_normal, directory = "GDCdata_TCGA_HNSC_RNASEQ_solid_tissue_normal")


# Data Preparation
# Extract and format the data for analysis
extract_data <- function(data) {
  df <- data@assays@data@listData[["tpm_unstrand"]]
  colnames(df) <- data@colData@listData[["patient"]]
  rownames(df) <- data@rowRanges@elementMetadata@listData[["gene_name"]]
  df <- rbind(df,
              race = data@colData@listData[["race"]],
              gender = data@colData@listData[["gender"]],
              stage = data@colData@listData[["ajcc_clinical_stage"]],
              sample_type = data@colData@listData[["sample_type"]])
  as.data.frame(df)
}

nt_df <- extract_data(hnsc_rnaseq_nt) # solid tissue normal
tp_df <- extract_data(hnsc_rnaseq_tp) # primary tumours

# Combine and transpose the data
full_table <- cbind(nt_df, tp_df)
transposed_full_table <- as.data.frame(t(full_table))
transposed_full_table <- na.omit(transposed_full_table)


# Split data by race
df_white <- transposed_full_table[transposed_full_table$race == "white", ]
df_black <- transposed_full_table[transposed_full_table$race == "black or african american", ]
df_asian <- transposed_full_table[transposed_full_table$race == "asian", ]

# Further split by sample type
df_white_Normal <- df_white[df_white$sample_type == "Solid Tissue Normal",]
df_white_tumor <- df_white[df_white$sample_type == "Primary Tumor",]

df_black_Normal <- df_black[df_black$sample_type == "Solid Tissue Normal",]
df_black_tumor <- df_black[df_black$sample_type == "Primary Tumor",]

df_asian_Normal <- df_asian[df_asian$sample_type == "Solid Tissue Normal",]
df_asian_tumor <- df_asian[df_asian$sample_type == "Primary Tumor",]

# Define stages and counts
stage <- c('Stage I', 'Stage II', 'Stage III', 'Stage IVA', 'Stage IVB', 'Stage IVC')
counts_asians <- c(0, 2, 2, 6, 0, 0)
counts_blacks <- c(1, 4, 6, 28, 5, 4)

# Repeat stages according to counts
samples_blacks <- rep(stage, counts_blacks)
samples_asians <- rep(stage, counts_asians)

# Select random samples for blacks and asians from white tumors
selected_samples_from_blacks <- do.call(rbind, lapply(samples_blacks, function(name) {
  matching_rows <- df_white_tumor[df_white_tumor$stage == name,]
  matching_rows[sample(nrow(matching_rows), 1), ]
}))

selected_samples_from_asians <- do.call(rbind, lapply(samples_asians, function(name) {
  matching_rows <- df_white_tumor[df_white_tumor$stage == name,]
  matching_rows[sample(nrow(matching_rows), 1), ]
}))

# Combine selected samples with original tumor data
df_w_a <- rbind(df_asian_tumor, selected_samples_from_asians)
df_w_b <- rbind(df_black_tumor, selected_samples_from_blacks)

# Optionally, if you need a sample from white Normal data
random_rows <- df_white_Normal[sample(nrow(df_white_Normal), 40), ]
random_rows$race <- paste0("white healthy")
df_w_w_b <- rbind(df_black_tumor, selected_samples_from_blacks, random_rows)

# Optionally, if you need a sample from white Normal data
random_rows <- df_white_Normal[sample(nrow(df_white_Normal), 10), ]
random_rows$race <- paste0("white healthy")
df_w_w_a <- rbind(df_asian_tumor, selected_samples_from_asians, random_rows)

# Now df_w_a and df_w_b contain the combined data for Asians and Blacks, respectively

all_ft <- apply(as.matrix.noquote(df_w_w_b[genes_w_w_b_bor$selected_features]), 2, as.numeric)
#all_ft <- apply(as.matrix.noquote(df_w_w_b[union_genes_w_b]), 2, as.numeric)
#all_ft <- apply(as.matrix.noquote(transposed_full_table[union_genes_w_b]), 2, as.numeric)
#all_ft <- scale(all_ft , center = T, scale = T)

all_ft <- log10(all_ft)
all_ft[is.infinite(all_ft)] <- -max(all_ft)
t_all_ft <- t(all_ft)
colnames(t_all_ft) <- paste("Sample", 1:ncol(t_all_ft), sep = "_")

race <- data.frame(df_w_w_b$race)
#stage <- data.frame(df_w_b$stage)

row.names(race) <- colnames(t_all_ft)
#row.names(stage) <- colnames(t_all_ft)

# Heatmap

heat_plot <- pheatmap(t_all_ft,
                      col = brewer.pal(10, 'RdYlGn'), # choose a colour scale for your data
                      cluster_rows = T, cluster_cols = T, # set to FALSE if you want to remove the dendograms
                      clustering_distance_cols = 'euclidean',
                      clustering_distance_rows = 'euclidean',
                      clustering_method = 'ward.D2',
                      cutree_rows = 9,
                      cutree_cols = 3,
                      #annotation_row = gene_functions_df, # row (gene) annotations
                      annotation_col = race, # column (sample) annotations
                      #annotation_colors = ann_colors, # colours for your annotations
                      #annotation_names_row = F, 
                      #annotation_names_col = F,
                      fontsize_row = 10,          # row label font size
                      fontsize_col = 7,          # column label font size 
                      angle_col = 45, # sample names at an angle
                      #legend_breaks = c(0, 5, 10), # legend customisation
                      #legend_labels = c("Low", "Medium", "High"), # legend customisation
                      show_colnames = F, show_rownames = F, # displaying column and row names
                      main = "Heatmap for Boruta Feature Selection Log transformed" # a title for our heatmap
                      ) 

########################################################################################

############################### Supervised Learning ####################################

########################################################################################
########################################################################################

##################################### Boruta ###########################################

########################################################################################

# Define a function to perform Boruta
feature_selection_boruta <- function(data, target) {
  set.seed(42)  # for reproducibility
  features <- apply(as.matrix.noquote(data), 2, as.numeric)
  
  target_encoded <- LabelEncoder$new()$fit_transform(target)
  boruta <- Boruta(features,
                   target_encoded,
                   pValue = 0.05,
                   maxRuns = 1000,
                   doTrace = 1)
  #plotImpHistory(
  #  boruta,
  #  colCode = c("green", "yellow", "red", "blue"),
  #  col = NULL,
  #  type = "l",
  #  lty = 1,
  #  pch = 0,
  #  xlab = "Classifier run",
  #  ylab = "Importance",
  #)
  #bor_imp_plot <- recordPlot()
  #dev.off()
  return(list(
    selected_features = getSelectedAttributes(boruta)
    #imp_plot = bor_imp_plot
  ))
}

# Predictions:
genes_w_a_bor <- feature_selection_boruta(df_w_a[, 1:60660], df_w_a$race)
genes_w_b_bor <- feature_selection_boruta(df_w_b[, 1:60660], df_w_b$race)
genes_w_w_b_bor <- feature_selection_boruta(df_w_w_b[, 1:60660], df_w_w_b$race)

genes_w_b_bor$selected_features

########################################################################################

#################################### XGBoost ###########################################

########################################################################################

# Define a function to perform XGBoost
feature_selection_xgboost <- function(data, target, cv_folds = 5, tune_grid = NULL) {
  # Ensure the target variable is a factor for binary classification
  data <- apply(as.matrix.noquote(data), 2, as.numeric)
  target_encoded <- LabelEncoder$new()$fit_transform(target)
  target <- as.factor(target_encoded)
  
  # Split the data into training and testing sets
  set.seed(42)  # for reproducibility
  trainIndex <- createDataPartition(target, p = 0.8, list = FALSE, times = 1)
  trainData <- data[trainIndex, ]
  testData <- data[-trainIndex, ]
  trainLabels <- target[trainIndex]
  testLabels <- target[-trainIndex]
  
  dtrain <- xgb.DMatrix(data = as.matrix(trainData), label = as.numeric(trainLabels) - 1)
  dtest <- xgb.DMatrix(data = as.matrix(testData), label = as.numeric(testLabels) - 1)
  
  # Define default parameter grid if none is provided
  if (is.null(tune_grid)) {
    tune_grid <- expand.grid(
      nrounds = c(100), 
      max_depth = c(3), 
      eta = c(0.1), 
      gamma = c(0, 1, 5), 
      colsample_bytree = c(1), 
      min_child_weight = c(3, 5),
      subsample = c(1)
    )
  }
  
  # Cross-validation
  train_control <- trainControl(method = "cv", number = cv_folds, verboseIter = TRUE, allowParallel = TRUE)
  xgb_tune <- train(x = trainData, y = trainLabels, method = "xgbTree", trControl = train_control, tuneGrid = tune_grid)
  best_params <- xgb_tune$bestTune
  
  # Train final model with the best parameters
  final_model <- xgb.train(params = list(
    objective = "binary:logistic", # logistic classification
    eval_metric = "logloss",
    max_depth = best_params$max_depth,
    eta = best_params$eta,
    gamma = best_params$gamma,
    colsample_bytree = best_params$colsample_bytree,
    min_child_weight = best_params$min_child_weight,
    subsample = best_params$subsample
  ), data = dtrain, nrounds = best_params$nrounds)
  
  # Predictions
  preds <- predict(final_model, dtest)
  preds <- ifelse(preds > 0.5, 1, 0)
  
  # Confusion matrix
  confusion <- confusionMatrix(as.factor(preds), testLabels)
  
  # Feature importance
  importance_matrix <- xgb.importance(model = final_model)
  xgb.plot.importance(importance_matrix)
  plot_imp_xgb <- recordPlot()
  dev.off()
  xgb_features <- importance_matrix$Feature
  
  return(list(
    best_params = best_params,
    final_model = final_model,
    confusion_matrix = confusion,
    feature_importance = importance_matrix,
    selected_features = xgb_features,
    plot_imp_xgb = plot_imp_xgb
  ))
}

# Predictions:
genes_w_a_xgb <- feature_selection_xgboost(data = df_w_a[, 1:60660], target = df_w_a$race) # not working properly due to small data size
genes_w_b_xgb <- feature_selection_xgboost(data = df_w_b[, 1:60660], target = df_w_b$race)
genes_w_w_b_xgb <- feature_selection_xgboost(data = df_w_w_b[, 1:60660], target = df_w_w_b$race)
genes_w_b_xgb$selected_features
########################################################################################

########################## Recursive Feature Elimination ###############################

########################################################################################

# Define a function to perform RFE
feature_selection_rfe <- function(data, target) {
  
  set.seed(42) # for reproducibility
  
  # Ensure the target variable is a factor for binary classification
  if (!is.factor(target)) {
    target <- as.factor(target)
  }
  
  # Convert data to numeric and remove constant columns
  data <- apply(as.matrix.noquote(data), 2, as.numeric)
  
  # Scale the data
  data <- scale(data)
  
  # Remove columns with NA values
  data <- data[, colSums(is.na(data)) == 0]
  
  # Set up the RFE control
  control <- rfeControl(functions = rfFuncs, # random forest
                        method = "repeatedcv", # repeated cv
                        repeats = 5, # number of repeats
                        number = 10,
                        verbose = TRUE) # number of folds
  
  # Run RFE with the tuned SVM model
  rfeProfile <- rfe(data,
                    target,
                    sizes = c(600),
                    rfeControl = control,
                    metric = 'Kappa'
                    )
  
  # Get the selected features
  rfe_features <- predictors(rfeProfile)
  
  varimp_data <- data.frame(feature = row.names(varImp(rfeProfile))[1:20], # top 20 only
                            importance = varImp(rfeProfile)[1:20, 1])
  
  rfe_imp_plot <- ggplot(data = varimp_data, 
         aes(x = reorder(feature, -importance), y = importance, fill = feature)) +
    geom_bar(stat="identity") + labs(x = "Features", y = "Variable Importance") + 
    geom_text(aes(label = round(importance, 2)), vjust=1.6, color="white", size=4) + 
    theme_bw() + theme(legend.position = "none")
  
  
  return(list(selected_features = rfe_features,
              rfe_imp_plot = rfe_imp_plot
              )
         )
}

# Predictions:
genes_w_a_rfe <- feature_selection_rfe(data = df_w_a[, 1:60660], target = df_w_a$race)
genes_w_b_rfe <- feature_selection_rfe(data = df_w_b[, 1:60660], target = df_w_b$race)
genes_w_w_b_rfe <- feature_selection_rfe(data = df_w_w_b[, 1:60660], target = df_w_w_b$race)

########################################################################################

################################## Random Forest #######################################

########################################################################################

# Define a function to perform Random Forest
feature_selection_rf <- function(data, target) {
  
  set.seed(42) # For reproducibility
  
  # Ensure the target variable is a factor for binary classification
  target<- as.factor(target)
  
  # Split the data into training and testing sets
  trainIndex <- createDataPartition(target, p = 0.8, list = FALSE)
  X_train <- data[trainIndex, ]
  y_train <- target[trainIndex]
  X_test <- data[-trainIndex, ]
  y_test <- target[-trainIndex]
  
  rf_model <- randomForest(x = X_train, y = y_train,
                           proximity=TRUE,
                           importance = TRUE)
  # Important features plot
  varImpPlot(rf_model)
  plot_imp_rf <- recordPlot()
  dev.off()
  # Extract the importance of each variable
  importance <- importance(rf_model)
  
  # Rank variables by importance
  importance_ranked <- importance[order(importance[, 1], decreasing = TRUE), ]
  
  # Select the most important features based on a threshold
  selected_features <-rownames(importance_ranked[importance_ranked[,'MeanDecreaseAccuracy'] != 0, ])
  
  
  
  # Optionally, use only the selected features for a new model
  X_train_selected <- X_train[, selected_features]
  X_test_selected <- X_test[, selected_features]
  
  # Fit a Random Forest model with selected features
  rf_model_selected <- randomForest(x = X_train_selected, y = y_train)
  
  # Evaluate the new model
  predictions <- predict(rf_model_selected, newdata = X_test_selected)
  
  # Classification metrics
  confusionMatrix <- confusionMatrix(predictions, y_test)
  
  # Plot ROC curve and calculate AUC
  prob_predictions <- predict(rf_model_selected, newdata = X_test_selected, type = "prob")
  roc_curve <- roc(y_test, prob_predictions[,2])
  
  plot(roc_curve)
  plot_roc <- recordPlot()
  dev.off()
  auc <- auc(roc_curve)
  
  
  return(list(
    plot_imp_rf = plot_imp_rf,
    importance_ranked = importance_ranked,
    selected_features = selected_features,
    confusion_matrix = confusionMatrix,
    roc = roc_curve,
    plot_roc = plot_roc,
    auc <- auc
  ))
}

# Predictions:
genes_w_a_rf <- feature_selection_rf(data = df_w_a[, 1:60660], target = df_w_a$race)
genes_w_b_rf <- feature_selection_rf(data = df_w_b[, 1:60660], target = df_w_b$race)
genes_w_w_b_rf <- feature_selection_rf(data = df_w_w_b[, 1:60660], target = df_w_w_b$race)
genes_w_b_rf$selected_features
########################################################################################

############################### Lasso Classification ###################################

########################################################################################

# Define a function to perform Lasso Classification
feature_selection_lasso <- function(data, target, alpha_value = 1) {
  
  data <- apply(as.matrix.noquote(data), 2, as.numeric)
  target <- LabelEncoder$new()$fit_transform(target)
  
  # Fit a Lasso regression model
  lasso_model <- glmnet(data, target, alpha = alpha_value, family = "multinomial") # binomial for classification
  
  # Perform cross-validation to find the optimal lambda
  cv_model <- cv.glmnet(data, target, alpha = alpha_value, family = "multinomial") # binomial for classification
  
  # Get the lambda that minimizes the cross-validation error
  optimal_lambda <- cv_model$lambda.min
  
  # Fit the Lasso model using the optimal lambda
  final_model <- glmnet(data, target, alpha = alpha_value, lambda = optimal_lambda)
  
  # Extract the coefficients of the final model
  coefficients <- coef(final_model)
  
  # Select the features with non-zero coefficients
  coefficients_matrix <- as.matrix(coefficients)
  selected_features <- rownames(coefficients_matrix)[coefficients_matrix != 0]
  lasso_features <- selected_features[selected_features != "(Intercept)"]  # Exclude the intercept
  
  return(list(optimal_lambda = optimal_lambda,
              coefficients = coefficients,
              selected_features = lasso_features))
}

# Predictions:
genes_w_a_lasso <- feature_selection_lasso(data = df_w_a[, 1:60660], target = df_w_a$race)
genes_w_b_lasso <- feature_selection_lasso(data = df_w_b[, 1:60660], target = df_w_b$race)
genes_w_w_b_lasso <- feature_selection_lasso(data = df_w_w_b[, 1:60660], target = df_w_w_b$race)
genes_w_b_lasso$selected_features

########################################################################################

#################################### Elastic Net #######################################

########################################################################################

# Define the function for Elastic Net
feature_selection_enet <- function(data, target, alpha_values = seq(0, 1, by = 0.1)) {
  
  data <- apply(as.matrix.noquote(data), 2, as.numeric)
  
  set.seed(42)

  # Split the data into training and testing sets
  trainIndex <- createDataPartition(target, p = 0.8, list = FALSE)
  X_train <- data[trainIndex, ]
  y_train <- target[trainIndex]
  X_test <- data[-trainIndex, ]
  y_test <- target[-trainIndex]
  
  # Perform cross-validation for each alpha value
  cv_results <- lapply(alpha_values, function(a) {
    cv.glmnet(X_train, y_train, alpha = a, family = "binomial")
  })
  
  # Find the best model based on cross-validated mean squared error
  best_model <- cv_results[[which.min(sapply(cv_results, function(mod) min(mod$cvm)))]]
  
  # Extract the best lambda and alpha
  best_lambda <- best_model$lambda.min
  best_alpha <- alpha_values[which.min(sapply(cv_results, function(mod) min(mod$cvm)))]
  
  # Train the final Elastic Net model with the best alpha and lambda
  final_model <- glmnet(X_train, y_train, alpha = best_alpha, lambda = best_lambda, family = "binomial")
  
  # Get the coefficients from the model
  coefficients <- coef(final_model)
  
  # Extract feature names with non-zero coefficients
  enet_features <- rownames(coefficients)[which(coefficients != 0)]
  enet_features <- enet_features[enet_features != "(Intercept)"]
  
  # Reduce the dataset to the selected features
  X_train_reduced <- X_train[, enet_features, drop = FALSE]
  X_test_reduced <- X_test[, enet_features, drop = FALSE]
  
  # Train a new Elastic Net model using the reduced feature set
  final_model_reduced <- glmnet(X_train_reduced, y_train, alpha = best_alpha, lambda = best_lambda, family = "binomial")
  
  # Predict and evaluate the model
  y_pred_train <- predict(final_model_reduced, X_train_reduced, type = "class")
  y_pred_test <- predict(final_model_reduced, X_test_reduced, type = "class")
  
  # Calculate evaluation metrics
  train_confusion <- confusionMatrix(as.factor(y_pred_train), as.factor(y_train))
  test_confusion <- confusionMatrix(as.factor(y_pred_test), as.factor(y_test))
  
  train_accuracy <- train_confusion$overall['Accuracy']
  test_accuracy <- test_confusion$overall['Accuracy']
  
  # Return results as a list
  results <- list(
    best_alpha = best_alpha,
    best_lambda = best_lambda,
    selected_features = enet_features,
    train_accuracy = train_accuracy,
    test_accuracy = test_accuracy,
    train_confusion_matrix = train_confusion$table,
    test_confusion_matrix = test_confusion$table
  )
  
  return(results)
}

# Predictions:
genes_w_a_enet <- feature_selection_enet(data = df_w_a[, 1:60660], target = df_w_a$race) # probably the model wrong because of sample size
genes_w_b_enet <- feature_selection_enet(data = df_w_b[, 1:60660], target = df_w_b$race)
genes_w_w_b_enet <- feature_selection_enet(data = df_w_w_b[, 1:60660], target = df_w_w_b$race)

########################################################################################

############################## Unsupervised Learning ###################################

########################################################################################
########################################################################################

########################## Principal Component Analysis ################################

########################################################################################

# Define the function for PCA
feature_selection_pca <- function(data, variance_threshold = 0.80, loading_threshold = 0.01) {
  
  set.seed(42)
  
  # Ensure the data frame is numeric
  numeric_data <- apply(as.matrix.noquote(data), 2, as.numeric)
  
  # Perform PCA
  pca_result <- prcomp(numeric_data)
  
  # Plot explained variance
  plot_imp_pca = fviz_eig(pca_result, addlabels = TRUE)
  
  # View summary of PCA results
  pca_summary <- summary(pca_result)
  
  # Determine the number of principal components that explain the specified amount of variance
  explained_variance <- pca_summary$importance[2, ]
  cumulative_variance <- cumsum(explained_variance)
  num_components <- which(cumulative_variance >= variance_threshold)[1]
  
  # Get the selected principal components
  selected_components <- pca_result$x[, 1:num_components]
  
  # Get the loadings
  loadings <- pca_result$rotation
  
  # Determine which features (genes) are significant based on the loadings
  selected_features <- abs(loadings[, 1:num_components]) > loading_threshold
  selected_genes <- rownames(loadings)[rowSums(selected_features) > 0]
  
  # Return the results as a list
  list(
    pca_result = pca_result,
    pca_summary = pca_summary,
    plot_imp_pca = plot_imp_pca,
    num_components = num_components,
    selected_components = selected_components,
    selected_features = selected_genes
  )
}

# Predictions:
genes_w_a_pca <- feature_selection_pca(data = df_w_a[, 1:60660], variance_threshold = 0.80, loading_threshold = 0.01)
genes_w_b_pca <- feature_selection_pca(data = df_w_b[, 1:60660], variance_threshold = 0.80, loading_threshold = 0.01)
genes_w_w_b_pca <- feature_selection_pca(data = df_w_w_b[, 1:60660], variance_threshold = 0.80, loading_threshold = 0.01)
genes_w_b_pca$selected_features

########################################################################################

########################### Pathway Enrichment Analysis ################################

########################################################################################

union_genes_w_a <- Reduce(union,
                          list(
                            genes_w_a_bor$selected_features,
                            #genes_w_a_xgb$selected_features, # not working due to sample size
                            genes_w_a_rfe$selected_features,
                            #genes_w_a_rf$selected_features, # to discuss with eleanor probably overfitting
                            genes_w_a_lasso$selected_features,
                            genes_w_a_enet$selected_features,
                            genes_w_a_pca$selected_features
                          )
)




union_genes_w_b <- Reduce(union,
                      list(
                        genes_w_b_bor$selected_features,
                        genes_w_b_xgb$selected_features,
                        genes_w_b_rfe$selected_features,
                        #genes_w_b_rf$selected_features, # to discuss with eleanor probably overfitting
                        genes_w_b_lasso$selected_features,
                        genes_w_b_enet$selected_features,
                        genes_w_b_pca$selected_features
                     )
)



union_genes_w_w_b <- Reduce(union,
                          list(
                            genes_w_w_b_bor$selected_features,
                            #genes_w_w_b_xgb$selected_features,
                            genes_w_w_b_rfe$selected_features,
                            #genes_w_w_b_rf$selected_features, # to discuss with eleanor probably overfitting
                            genes_w_w_b_lasso$selected_features,
                            #genes_w_b_enet$selected_features,
                            genes_w_w_b_pca$selected_features
                          )
)


#print(msigdbr_collections(), n=23)
#msigdbr_species()
#msigdbr(species = "Homo sapiens")

cgp_gene_sets = msigdbr(species = "Homo sapiens",
                        category = "C7",
                        subcategory = 'IMMUNESIGDB' 
                        )
msigdbr_t2g = cgp_gene_sets %>% dplyr::distinct(gs_name, gene_symbol) %>% as.data.frame()

enrichment <- enricher(gene = union_genes_w_a, TERM2GENE = msigdbr_t2g, pvalueCutoff = 0.01)


dotplot(enrichment, 
        color = "p.adjust", 
        font.size = 12, 
        showCategory=20,
        title = "MSigDB Pathway Enrichment oncogenic signature gene sets", 
        orderBy = "x"
        )


barplot(enrichment, 
        color = "p.adjust", 
        showCategory=20, 
        font.size = 12,
        title = "MSigDB Pathway Enrichment immunologic signature gene sets"
        )

cnetplot(enrichment,
         showCategory = 11,
         layout= 'kk',
         colorEdge = TRUE,
         node_label = 'category'
         )

cnetplot(enrichment, circular = TRUE, colorEdge = TRUE) 

mutate(enrichment, qscore = -log(p.adjust, base=10)) %>% 
  barplot(x="qscore",showCategory=10, font.size = 12) + ggtitle("MSigDB Pathway Enrichment oncogenic signature gene sets as Log 10 P adjusted values")

# Get the similarity matrix
enrichmentx <- pairwise_termsim(enrichment)

emapplot(enrichmentx, node_label = 'none') + ggtitle("MSigDB Enrichment Map")



# Similarity space plot of enrichment analysis results
ssplot(enrichmentx,
       showCategory = 100,
       repel = T,
       color = "p.adjust",
       node_label = 'none'
       )
treeplot(enrichmentx,
         showCategory = 10, 
         color = "p.adjust", 
         cluster.params = list(method ="average"),
         )

upsetplot(enrichment)



# Define the function
perform_enrichment_analysis <- function(selected_features, msigdbr_t2g, pvalue_cutoff = 0.01) {
  
  # Perform enrichment analysis
  enrichment <- enricher(gene = selected_features, TERM2GENE = msigdbr_t2g, pvalueCutoff = pvalue_cutoff)
  
  # Initialize an empty list to store plots
  plot_list <- list()
  
  # Dotplot
  dotplot_plot <- dotplot(enrichment, 
                          color = "p.adjust", 
                          font.size = 10, 
                          title = "MSigDB Pathway Enrichment immunologic signature gene sets", 
                          orderBy = "x")
  plot_list$dotplot <- dotplot_plot
  
  # Barplot
  barplot_plot <- barplot(enrichment, 
                          color = "p.adjust", 
                          showCategory = 10, 
                          font.size = 8,
                          title = "MSigDB Pathway Enrichment immunologic signature gene sets")
  plot_list$barplot <- barplot_plot
  
  # Cnetplot with layout 'kk'
  cnetplot_plot <- cnetplot(enrichment,
                            showCategory = 4,
                            layout = 'kk')
  plot_list$cnetplot <- cnetplot_plot
  
  # Calculate qscore and plot
  barplot_qscore <- mutate(enrichment, qscore = -log(p.adjust, base=10)) %>% 
    barplot(x="qscore", font.size = 12)
  plot_list$barplot_qscore <- barplot_qscore
  
  # Similarity matrix
  enrichmentx <- pairwise_termsim(enrichment)
  
  # Similarity space plot
  ssplot_plot <- ssplot(enrichmentx,
                        showCategory = 30)
  plot_list$ssplot <- ssplot_plot
  
  # Treeplot
  treeplot_plot <- treeplot(enrichmentx,
                            showCategory = 10, 
                            color = "p.adjust", 
                            cluster.params = list(method = "average"))
  plot_list$treeplot <- treeplot_plot
  
  # Upset plot
  upsetplot_plot <- upsetplot(enrichment)
  plot_list$upsetplot <- upsetplot_plot
  
  return(plot_list)
}

cgp_gene_sets = msigdbr(species = "Homo sapiens",
                        category = "C6"
)
msigdbr_t2g <- cgp_gene_sets %>% dplyr::distinct(gs_name, gene_symbol) %>% as.data.frame()
plots <- perform_enrichment_analysis(selected_features = union_genes_w_b, msigdbr_t2g = msigdbr_t2g)









##### MAF


# Somatic mutation data
query_mutation <- GDCquery(
  project = "TCGA-HNSC", 
  data.category = "Simple Nucleotide Variation",
  data.type = "Masked Somatic Mutation", 
  workflow.type = "Aliquot Ensemble Somatic Variant Merging and Masking"
)

View(getResults(query_mutation))

GDCdownload(
  query_mutation,
  token.file,
  method = "api",
  directory = "GDCdata_TCGA_HNSC_mutation",
  files.per.chunk = NULL
)

maf <- GDCprepare(query_mutation, directory = "GDCdata_TCGA_HNSC_mutation")

maf <- read.maf(maf = maf)

#Shows sample summry.
getSampleSummary(maf)
#Shows gene summary.
getGeneSummary(maf)
#shows clinical data associated with samples
getClinicalData(maf)
#Shows all fields in MAF
getFields(maf)

# Plotting MAF summary
plotmafSummary(maf = asians_maf, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)

#oncoplot for top ten mutated genes
oncoplot(maf = whites_maf, top = 10)


# Transition and Transversions.
hnsc.titv = titv(maf = maf, plot = FALSE, useSyn = TRUE)
#plot titv summary
plotTiTv(res = hnsc.titv)

# Compare mutation load against TCGA cohorts
hnsc.mutload = tcgaCompare(maf = whites_maf, cohortName = 'HNSC', logscale = TRUE, capture_size = 50)

#Plotting Variant Allele Frequencies VAF
plotVaf(maf = maf, vafCol = 'i_TumorVAF_WU')


# Somatic Interactions (which performs pair-wise Fisher’s Exact test to detect such significant pair of genes)
# exclusive/co-occurance event analysis on top 10 mutated genes.
somaticInteractions(maf = asians_maf, 
                    top = 40, 
                    pvalue = c(0.05, 0.01),
                    fontSize = 0.8,
                    leftMar = 8,
                    topMar = 7,
                    #showCounts = T,
                    sigSymbolsSize = 5,
                    sigSymbolsFontSize = 1.5,
                    #pvSymbols = c(46, 42)
                    )


# Detecting cancer driver genes based on positional clustering
hnsc.sig = oncodrive(maf = blacks_maf, AACol = 'Protein_Change', minMut = 5, pvalMethod = 'zscore')
plotOncodrive(res = hnsc.sig, fdrCutOff = 0.1, useFraction = TRUE, labelSize = 0.5)




blacks_index <- which(hnsc_rnaseq_tp@colData@listData[["race"]] == 'black or african american')
whites_index <- which(hnsc_rnaseq_tp@colData@listData[["race"]] == 'white')
asians_index <- which(hnsc_rnaseq_tp@colData@listData[["race"]] == 'asian')

blacks_code_bar <- as.list(hnsc_rnaseq_tp@colData@rownames[blacks_index])
whites_code_bar <- as.list(hnsc_rnaseq_tp@colData@rownames[whites_index])
asians_code_bar <- as.list(hnsc_rnaseq_tp@colData@rownames[asians_index])

# Function to extract the specific part of the ID
extract_id_part <- function(id) {
  sub("^([^-]+-[^-]+-[^-]+).*", "\\1", id)
}

# Extract the required parts from both lists
extracted_list1 <- sapply(blacks_code_bar, extract_id_part)
extracted_list2 <- sapply(whites_code_bar, extract_id_part)
extracted_list3 <- sapply(asians_code_bar, extract_id_part)

# Extract elements from list2 that have the same prefix as any element in extracted_list1
matching_elements <- maf@clinical.data[["Tumor_Sample_Barcode"]][sapply(maf@clinical.data[["Tumor_Sample_Barcode"]], function(id2) {
  extracted_part2 <- extract_id_part(id2)
  extracted_part2 %in% extracted_list1
})]

blacks_maf <- subsetMaf(maf = maf, 
                        tsb = matching_elements,
                        genes = union_genes_w_b
                        )

matching_elements2 <- maf@clinical.data[["Tumor_Sample_Barcode"]][sapply(maf@clinical.data[["Tumor_Sample_Barcode"]], function(id2) {
  extracted_part2 <- extract_id_part(id2)
  extracted_part2 %in% extracted_list2
})]
####
resultttt <- substr(gsub("\\.", "-", rownames(selected_samples_from_blacks)), 1, 12)
matching_elements23 <- maf@clinical.data[["Tumor_Sample_Barcode"]][sapply(maf@clinical.data[["Tumor_Sample_Barcode"]], function(id2) {
  extracted_part2 <- extract_id_part(id2)
  extracted_part2 %in% resultttt
})]
####


whites_maf <- subsetMaf(maf = maf, 
                        tsb = matching_elements23,
                        genes = union_genes_w_b
                        )

matching_elements3 <- maf@clinical.data[["Tumor_Sample_Barcode"]][sapply(maf@clinical.data[["Tumor_Sample_Barcode"]], function(id2) {
  extracted_part2 <- extract_id_part(id2)
  extracted_part2 %in% extracted_list3
})]

asians_maf <- subsetMaf(maf = maf, 
                        tsb = matching_elements3,
                        genes = union_genes_w_b)



hnc_b_w <- mafCompare(m1 = blacks_maf, m2 = whites_maf, m1Name = 'Blacks', m2Name = 'Whites', minMut = 5)
forestPlot(mafCompareRes = hnc_b_w, pVal = 0.05)

#whites top features
df_w_w_b[c("COL17A1", "TMCC1", 'LAMB3', 'FN1', 'FABP4', 'DSP', 'COL1A2', 'CD74', 'CD209', 'BGN')]
#blacks top features
df_w_w_b[c("LAMC2", "COL28A1", 'CD209', 'AKR1B10', 'DSP', 'COL3A1', 'COL1A2', 'CDK11A', 'CCDC144A', 'B2M')]


df_w_w_a[c('KRAS', 'HRAS', 'NRAS', 'EIF4G1', 'RB1', 'SHH', 'P53')]

df_w_w_a[c('XKR3', 'SFN', 'KRT5', 'JUP', 'DSP', 'COL1A2', 'COL17A1')]

df_w_w_b[c('KRAS', 'HRAS', 'NRAS', 'EIF4G1', 'RB1', 'SHH', 'ERBB2')]

features_select <- df_ordered_variance
features_select <- cbind(features_select, df_w_w_a$race)

ft <- apply(as.matrix.noquote(features_select[, -ncol(features_select)]), 2, as.numeric)

subset(features_select, select = 'df_w_w_a$race')

features_select <- features_select %>%
  rename(
    race = 'df_w_w_a$race'
  )

featurePlot(x = ft, 
            y = as.factor(features_select$race), 
            plot = "density", 
            ## Pass in options to xyplot() to 
            ## make it prettier
            scales = list(x = list(relation="free"), 
                          y = list(relation="free")), 
            adjust = 1.5, 
            pch = "|", 
            layout = c(4, 1), 
            auto.key = list(columns = 3))


featurePlot(x = ft[,1:10], 
            y = as.factor(features_select$race), 
            plot = "box", 
            ## Pass in options to bwplot() 
            scales = list(y = list(relation="free"),
                          x = list(rot = 90)),  
            layout = c(5,2 ), 
            auto.key = list(columns = 2))




df_w_w_b_union <- df_w_w_b[union_genes_w_b]
df_w_w_a_union <- df_w_w_a[union_genes_w_a]

df_w_b_union <- df_w_b[union_genes_w_b]
df_w_a_union <- df_w_a[union_genes_w_a]
# Calculate variances
variances <- apply(df_w_a_union, 2, var)

# Get ordered indices of columns based on variances
ordered_indices <- order(variances, decreasing = TRUE)

# Reorder data frame by variance
df_ordered_variance <- df_w_w_a_union[, ordered_indices]


mf <- apply(as.matrix.noquote(df_w_b_union), 2, as.numeric)

# Calculate means of each column
means <- colMeans(mf)

# Get ordered indices of columns based on means
ordered_indices <- order(means, decreasing = TRUE)

# Reorder data frame by mean
df_ordered_mean <- df_w_w_b_union[, ordered_indices]



# Calculate standard deviation for each column
std_devs <- apply(df_w_b_union, 2, sd)

# Get ordered indices of columns based on standard deviations
ordered_indices <- order(std_devs, decreasing = TRUE)

# Reorder data frame by standard deviation
df_ordered_std <- df_w_b_union[, ordered_indices]







df_w_b_union <-cbind(df_w_b_union, df_w_b$race)

