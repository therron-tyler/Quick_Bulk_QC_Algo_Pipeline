# SHAPES: Shape #'s 21 through 25 can be filled using specified colors within the ggplot function
# CONDITIONS FILE: 4 columns, No headers. 

# Tyler Therron 
# Winter Lab - Macrophage Genomics
# Northwestern University - Feinberg School of Medicine, Department of Rheumatology

# PCA_script_v3_interactive.R
# USAGE: Rscript PCA_script_v3_interactive.R data_filename.txt conditions_filename.csv "Experimental Groups"

libs <- .libPaths("/home/ttm3567/63_tylert/Analysis_Algorithms/Histogram_PCA_Renv/")


# Load necessary libraries
library(ggplot2)
library(dplyr)
library(rlang)
library(plotly)
library(htmlwidgets)

## Read command line arguments
args <- commandArgs(TRUE)
file_name <- args[1]
conds_file <- read.csv(args[2], header=FALSE)
sample_groups <- as.character(args[3])
output_directory <- as.character(args[4])

# file_name <- "/Users/ttm3567/Documents/November2024/BMC_Timecourse/BMC_Timecourse_HQsamples_CPM.txt"
# conds_file <- read.csv("/Users/ttm3567/Documents/November2024/BMC_Timecourse/Cell_Type_analyses/neg_don_PCA_sample_list.csv", header = FALSE)
# sample_groups <- as.character("BMC_Full_Timecourse_NegDon_CellType")
# output_directory <- "/Users/ttm3567/Documents/November2024/BMC_Timecourse/Cell_Type_analyses/"

# Read counts data
data <- read.table(file_name, header=TRUE, sep="\t", row.names = 1, check.names = FALSE)
colnames(data)
colnames(conds_file)
data$Gene_Symbol <- NULL

# Assign column names
colnames(conds_file) <- c("SampleID", "Condition", "Color", "Shape")
conds_file
## Data cleaning and verification with samples in conditions file
namelist <- conds_file$SampleID
print(namelist)

print(colnames(data))

missing_names <- setdiff(namelist, colnames(data))

# Print the missing names
if(length(missing_names) > 0) {
  cat("The following names are in namelist but not in colnames(data):\n")
  print(missing_names)
} else {
  cat("All names in namelist are present in colnames(data).\n")
}

namelist <- namelist[namelist %in% colnames(data)]
nrow(conds_file)

conds_file <- conds_file[conds_file$SampleID %in% colnames(data), ]
nrow(conds_file)

data <- data[ ,match(namelist, colnames(data))]

head(data)

# Define group conditions, print user-defined experimental conditions, colors, and shapes to the console
conds <- conds_file$Condition
for(i in 1:length(conds)) {
	  print(paste0("Experimental conditions [", i, "/", nrow(conds_file), "] -> ", conds[i]))
}
color_mapping <- conds_file$Color
for(i in 1:length(color_mapping)) {
	  print(paste0("Colors [", i, "/", nrow(conds_file), "] -> ", color_mapping[i]))
}

shape_mapping <- conds_file$Shape
for(i in 1:length(shape_mapping)) {
	  print(paste0("Shapes [", i, "/", nrow(conds_file), "] -> ", shape_mapping[i]))
}



## Filter out genes with zero counts
data <- data.matrix(data)
rs <- rowSums(data)
use <- (rs > 0)
data <- data[use, ]

## Define file prefix
file_prefix <- sample_groups

## Transpose the data and convert to matrix
data_tr <- t(data)
rownames(data_tr) <- colnames(data)

### PCA calculation
pca_comp <- prcomp(data_tr, scale. = TRUE, center = TRUE)
percentVar <- pca_comp$sdev^2/sum(pca_comp$sdev^2)

## Create DataFrame with PCA loadings for plotting
# Get the conditions in the desired order
desired_condition_order <- unique(conds_file$Condition)

pca_df <- data.frame(
		       PC1 = pca_comp$x[,1],
		         PC2 = pca_comp$x[,2],
			   conds,
			     sampleLab = rownames(pca_comp$x),
			       shape_mapping,
			         check.names = FALSE
		       )
colnames(pca_df)[3] <- sample_groups
colnames(pca_df)[3]
print(pca_df)

pca_df$shape_mapping <- as.factor(pca_df$shape_mapping)

pca_df$combined_factor <- paste0(pca_df[[sample_groups]], "_", pca_df$shape_mapping)
pca_df$combined_factor <- as.factor(pca_df$combined_factor)


# Generate color and shape mappings
color_mapping_plot <- setNames(color_mapping, pca_df$combined_factor)
shape_mapping_plot <- setNames(shape_mapping, pca_df$combined_factor)


# Map ggplot2 shape codes to Plotly symbol names
shape_mapping_ggplot_to_plotly <- c(
  "21" = "circle",          # Hollow circle with fill
  "22" = "square",          # Hollow square with fill
  "23" = "diamond",         # Hollow diamond with fill
  "24" = "triangle-up",     # Hollow upward triangle with fill
  "25" = "triangle-down",   # Hollow downward triangle with fill
  "0"  = "square-open",     # Hollow square
  "1"  = "circle-open",     # Hollow circle
  "2"  = "triangle-up-open",# Hollow upward triangle
  "5"  = "diamond-open",    # Hollow diamond
  "6"  = "triangle-down-open",# Hollow downward triangle
  "15" = "square",          # Solid square
  "16" = "circle",          # Solid circle
  "17" = "triangle-up",     # Solid upward triangle
  "18" = "diamond"          # Solid diamond
)

# Create a new column in pca_df for Plotly symbols
pca_df$plotly_symbol <- shape_mapping_ggplot_to_plotly[as.character(pca_df$shape_mapping)]
pca_df

# Generate custom legend labels from experimental group column
condition_labels <- unique(pca_df[[sample_groups]])
combined_factors <- unique(pca_df$combined_factor)
label_mapping <- setNames(condition_labels, combined_factors)
# label_mapping <- sapply(label_mapping, function(x) gsub("_", " ", x))
label_mapping

pca_df$combined_factor <- gsub("_[0-9]+$", "", pca_df$combined_factor)
names(color_mapping_plot) <- gsub("_[0-9]+$", "", names(color_mapping_plot))
names(label_mapping) <- gsub("_[0-9]+$", "", names(label_mapping))
pca_df$combined_factor
names(color_mapping_plot)
names(label_mapping)
# pca_df$combined_factor

for(i in 1:length(label_mapping)) {
	  print(paste0("Label Mapping [", i, "/", length(label_mapping), "] -> ", label_mapping[i]))
}


## PCA plot using Plotly
# Create the ggplot object
p <- pca_df %>%
  ggplot(aes(
	         x = PC1,
		     y = PC2,
		         label = sampleLab,
			     fill = combined_factor,
			         shape = combined_factor,
				     text = paste(
						        "Sample:", sampleLab,
							      "<br>Condition:", pca_df[[sample_groups]],
							            "<br>PC1:", round(PC1, 2),
								          "<br>PC2:", round(PC2, 2)
							    )
		   )) +
  geom_point(size = 5) +
    scale_shape_manual(
		         values = as.numeric(shape_mapping_plot),
			       labels = label_mapping,
			       name = sample_groups
			     ) +
  scale_fill_manual(
		      values = color_mapping_plot,
			    labels = label_mapping,
			    name = sample_groups
			  ) +
  labs(
           x = paste0("PC1: ", round(percentVar[1] * 100, 1), "% Variance Explained"),
	       y = paste0("PC2: ", round(percentVar[2] * 100, 1), "% Variance Explained"),
	           fill = sample_groups,
		       shape = sample_groups
	     ) +
  theme_linedraw(base_size = 16) +
    theme(
	      panel.grid.major = element_blank(),
	          panel.grid.minor = element_blank(),
		      legend.title = element_text(size = 15),
		          legend.text = element_text(size = 12),
			      legend.position = "right",
			          axis.text = element_text(size = 12),
				      axis.title = element_text(size = 18)
	        ) +
  guides(color = guide_legend(ncol = 1))

  # Convert ggplot to plotly object for interactivity
  p_plotly <- ggplotly(p, tooltip = "text")
  
  # Update the marker symbols in the plotly object
  for (i in seq_along(p_plotly$x$data)) {
	    if (p_plotly$x$data[[i]]$type == "scatter") {
		        # Get the combined_factor for this trace
		        trace_name <- p_plotly$x$data[[i]]$name
		        trace_name
      # Get the corresponding Plotly symbol
      symbol <- unique(pca_df$plotly_symbol[pca_df$combined_factor == trace_name])
          # Update the marker symbol
          p_plotly$x$data[[i]]$marker$symbol <- symbol
        }
  }
  
  desired_order <- unique(conds_file$Condition)
  
  # Extract trace names and indices
  trace_names <- sapply(p_plotly$x$data, function(trace) trace$name)
  trace_indices <- seq_along(p_plotly$x$data)
  
  # Create a data frame mapping trace names to indices
  trace_df <- data.frame(
    index = trace_indices,
    name = trace_names,
    stringsAsFactors = FALSE
  )
  
  # Map desired order to trace indices
  desired_trace_indices <- sapply(desired_order, function(name) {
    idx <- which(trace_df$name == name)
    if (length(idx) == 0) {
      NA  # Handle conditions not present in the plot
    } else {
      idx
    }
  })
  
  # Remove NAs (if any conditions are not present in the plot)
  desired_trace_indices <- desired_trace_indices[!is.na(desired_trace_indices)]
  p_plotly$x$data <- p_plotly$x$data[desired_trace_indices]
  
  # Extract the names again to verify
  new_trace_names <- sapply(p_plotly$x$data, function(trace) trace$name)
  print(new_trace_names)
  
  p_plotly
  # Save the interactive plot as an HTML file
  html_file_name <- paste0(output_directory, file_prefix, "_PCAplot.html")
  saveWidget(p_plotly, file = html_file_name, selfcontained = TRUE)

  ## Write CSV file with PCA data
  write.csv(pca_df, paste0(output_directory, file_prefix, "_PCA_Data.csv"), row.names = FALSE)

  # Print completion message
  print(paste0("Interactive PCA plot saved as ", html_file_name))

  
