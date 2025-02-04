# Tyler Therron 
# Winter Lab - Macrophage Genomics
# Northwestern University - Feinberg School of Medicine, Department of Rheumatology


# Usage: Rscript Histogram_script.R CPM_SLE1_3_MIC_trim.txt conditions.csv
# CONDITIONS FILE: 4 columns, No headers. 

libs <- .libPaths("/home/ttm3567/63_tylert/Analysis_Algorithms/Histogram_PCA_Renv/")
#install.packages("readr")

library(reshape2) 
library(ggplot2)
library(readr)
library(tidyr)
library(plotly)
library(htmlwidgets)
library(dplyr)

args <- commandArgs(TRUE)

data <- read_tsv(args[1])
conds_file <- read_csv(args[2], col_names=FALSE)
name <-  as.character(args[3])
output_directory <- as.character(args[4])

# Functions ----------- # ----------- # ------------
assign_names_and_reorder <- function(dataframe, conds_file) {
  print("assigning Ensembl IDs to be rownames")
  rownames(data) <- unlist(dataframe[,1])
  data[,1] <- NULL
  data <- data %>% select(-Gene_Symbol)
  
  print("reorder columns")
  data <- data[, conds_file$X1]
  
  print("filter by conditions file")
  data <- data[, colnames(data) %in% conds_file$X1]
  
  return(data)
}

log_transform_df_clean <- function(dataframe) {
  # Identify the numeric columns
  print("converting data to a log2 scale and adding 1")
  df <- dataframe
  numeric_cols <- sapply(df, is.numeric)
  df[, numeric_cols] <- log2(df[, numeric_cols] + 1)
  
  print("filter out rows with no gene expression")
  df <- df[rowSums(df) != 0, ]
  return(df)
}

melt_and_merge <- function(df) {
  # Melt the data frame
  df_melt <- melt(df, variable.name = "sample", value.name = "expression")
  
  # Merge with conditions
  df_melt <- merge(df_melt, conditions_ordered, by = "sample")
  return(df_melt)
}


# Functions ----------- # ----------- # ------------
# Predefined Palette ----------- # ----------- # ------------
color_dev = c("#000000", "#FFFF00", "#1CE6FF", "#FF34FF", "#FF4A46","#008941", 
                       "#006FA6", "#A30059", "#7A4900", "#0000A6", "#63FFAC", "#B79762",
                       "#004D43", "#8FB0FF", "#997D87","#5A0007", "#809693", "#1B4400", 
                       "#4FC601", "#3B5DFF", "#4A3B53", "#FF2F80","#61615A", "#BA0900", 
                       "#6B7900", "#00C2A0", "#FFAA92", "#FF90C9", "#B903AA", "#D16100",
                       "#000035", "#7B4F4B", "#A1C299", "#300018", "#0AA6D8", "#013349", 
                       "#00846F","#372101", "#FFB500", "#C2FFED", "#A079BF", "#CC0744", 
                       "#C0B9B2", "#001E09","#00489C", "#6F0062", "#0CBD66", "#EEC3FF", 
                       "#456D75", "#B77B68", "#7A87A1", "#788D66","#885578", "#FAD09F", 
                       "#FF8A9A", "#D157A0", "#BEC459", "#456648", "#0086ED", "#886F4C",
                       "#34362D", "#B4A8BD", "#00A6AA", "#452C2C", "#636375", "#A3C8C9", 
                       "#FF913F", "#938A81","#575329", "#00FECF", "#B05B6F", "#8CD0FF", 
                       "#3B9700", "#04F757", "#C8A1A1", "#1E6E00","#7900D7", "#A77500", 
                       "#6367A9", "#A05837", "#6B002C", "#772600", "#D790FF", "#9B9700",
                       "#549E79", "#FFF69F", "#201625", "#72418F", "#BC23FF", "#99ADC0", 
                       "#3A2465", "#922329","#5B4534", "#FDE8DC", "#404E55", "#0089A3", 
                       "#CB7E98", "#A4E804", "#324E72", "#6A3A4C","#83AB58", "#001C1E", 
                       "#D1F7CE", "#004B28", "#C8D0F6", "#A3A489", "#806C66", "#222800",
                       "#BF5650", "#E83000", "#66796D", "#DA007C", "#FF1A59", "#8ADBB4", 
                       "#1E0200", "#5B4E51","#C895C5", "#320033", "#FF6832", "#66E1D3", 
                       "#CFCDAC", "#D0AC94", "#7ED379", "#012C58","#7A7BFF", "#D68E01", 
                       "#353339", "#78AFA1", "#FEB2C6", "#75797C", "#837393", "#943A4D",
                       "#B5F4FF", "#D2DCD5", "#9556BD", "#6A714A", "#001325", "#02525F", 
                       "#0AA3F7", "#E98176","#DBD5DD", "#5EBCD1", "#3D4F44", "#7E6405", 
                       "#02684E", "#962B75", "#8D8546", "#9695C5","#E773CE", "#D86A78", 
                       "#3E89BE", "#CA834E", "#518A87", "#5B113C", "#55813B", "#E704C4",
                       "#00005F", "#A97399", "#4B8160", "#59738A", "#FF5DA7", "#F7C9BF", 
                       "#643127", "#513A01","#6B94AA", "#51A058", "#A45B02", "#1D1702", 
                       "#E20027", "#E7AB63", "#4C6001", "#9C6966","#64547B", "#97979E", 
                       "#006A66", "#391406", "#F4D749", "#0045D2", "#006C31", "#DDB6D0",
                       "#7C6571", "#9FB2A4", "#00D891", "#15A08A", "#BC65E9", "#FFFFFE", 
                       "#C6DC99", "#203B3C","#671190", "#6B3A64", "#F5E1FF", "#FFA0F2", 
                       "#CCAA35", "#374527", "#8BB400", "#797868","#C6005A", "#3B000A", 
                       "#C86240", "#29607C", "#402334", "#7D5A44", "#CCB87C", "#B88183",
                       "#AA5199", "#B5D6C3", "#A38469", "#9F94F0", "#A74571", "#B894A6", 
                       "#71BB8C", "#00B433","#789EC9", "#6D80BA", "#953F00", "#5EFF03", 
                       "#E4FFFC", "#1BE177", "#BCB1E5", "#76912F","#003109", "#0060CD", 
                       "#D20096", "#895563", "#29201D", "#5B3213", "#A76F42", "#89412E",
                       "#1A3A2A", "#494B5A", "#A88C85", "#F4ABAA", "#A3F3AB", "#00C6C8", 
                       "#EA8B66", "#958A9F","#BDC9D2", "#9FA064", "#BE4700", "#658188", 
                       "#83A485", "#453C23", "#47675D", "#3A3F00","#061203", "#DFFB71", 
                       "#868E7E", "#98D058", "#6C8F7D", "#D7BFC2", "#3C3E6E", "#D83D66",
                       "#2F5D9B", "#6C5E46", "#D25B88", "#5B656C", "#00B57F", "#545C46", 
                       "#866097", "#365D25","#252F99", "#00CCFF", "#674E60", "#FC009C","#FFDBE5", "#92896B")
# Predefined Palette ----------- # ----------- # ------------                      

# data <- read_tsv("/Users/ttm3567/Documents/November2024/BMC_Timecourse/BMC_Timecourse_HQsamples_CPM.txt")
# conds_file <- read_csv("/Users/ttm3567/Documents/November2024/BMC_Timecourse/BMC_Timecourse_PCA_sample_list.csv", col_names = FALSE)
# name <- "BMC_Full_Timecourse"
# output_directory <- "/Users/ttm3567/Documents/November2024/BMC_Timecourse/"

print("Data: ")
print(data[1:10,])

print("Conditions: ")
print(conds_file)

nrow(conds_file)
conds_file <- conds_file[conds_file$X1 %in% colnames(data), ]
nrow(conds_file)

print("converted to dataframe after read in of data and conditions")
data <- as.data.frame(data)

data <- assign_names_and_reorder(data, conds_file)

df <- log_transform_df_clean(data)

# Create a Data Frame Mapping Samples to Conditions
conditions_ordered <- data.frame(
  sample = conds_file$X1,
  conds = conds_file$X2,
  stringsAsFactors = FALSE
)

# Assign a Numeric Identifier to Each Unique Condition
conditions_ordered$cond_id <- as.numeric(factor(conditions_ordered$conds))

# Create a Mapping from Condition Names to IDs
condition_mapping <- unique(conditions_ordered[, c("conds", "cond_id")])

print("Experimental Group to ID mapping")
print(condition_mapping)
print("Conditions and their associated samples:")
print(split(conditions_ordered$sample, conditions_ordered$conds))

df_melt <- melt_and_merge(df)

# Define bin edges
bin_edges <- seq(0, max(df_melt$expression) + 0.5, by = 0.5)

# Cut the expression data into bins
df_melt$bin <- cut(df_melt$expression, breaks = bin_edges, include.lowest = TRUE, right = FALSE)

# Count occurrences in each bin for each sample
bins_pt5 <- df_melt %>%
  group_by(sample, bin) %>%
  summarize(count = n()) %>%
  ungroup()

# Create a wide-format data frame with bins as rows and samples as columns
bins_wide <- bins_pt5 %>%
  pivot_wider(names_from = sample, values_from = count, values_fill = 0)

# Convert to a matrix for easier manipulation
bins_matrix <- as.matrix(bins_wide[,-1])
rownames(bins_matrix) <- bins_wide$bin

# Calculate cumulative counts for each sample
re_bins_pt5 <- apply(bins_matrix, 2, function(x) {
  rev(cumsum(rev(x)))
})

# Convert back to data frame
re_bins_pt5 <- as.data.frame(re_bins_pt5)
rownames(re_bins_pt5) <- bins_wide$bin
bins_wide$bin

re_bins_pt5
# Extract the lower limit of each bin and format the row names
rownames(re_bins_pt5) <- sapply(strsplit(rownames(re_bins_pt5), ","), function(x) {
  paste0(substr(x[1], 2, nchar(x[1])), " <")
})

# Write binned counts to CSV
write.csv(bins_matrix, file = paste0(output_directory, name, "_bincounts_0.5.csv"), row.names = TRUE)

# Write cumulative counts to CSV
write.csv(re_bins_pt5, file = paste0(output_directory, name, "_incremental_bincounts_0.5.csv"), row.names = TRUE)

# Prepare data for plotting
re_bins_pt5$bin <- rownames(re_bins_pt5)
df_plot <- re_bins_pt5 %>%
  pivot_longer(-bin, names_to = "sample", values_to = "cumulative_count")

# Convert 'bin' to a numeric value representing the lower edge of each bin
df_plot$bin_numeric <- as.numeric(gsub(" <", "", df_plot$bin))

# Merge with conditions to get sample colors
df_plot <- merge(df_plot, conditions_ordered[, c("sample", "conds")], by = "sample")

# Generate color palette
num_samples <- length(unique(df_plot$sample))
num_samples
palette_colors <- color_dev[1:num_samples]
palette_colors

# Create the plot
cum_count <- ggplot(df_plot, aes(x = bin_numeric, y = cumulative_count, color = sample)) +
  geom_line() +
  scale_color_manual(values = palette_colors) +
  theme_minimal() +
  labs(
    title = "Histogram of Gene Expression Counts",
    x = "Bins (Size = 0.5)",
    y = "Cumulative Count",
    color = "Sample"
  ) +
  coord_cartesian(ylim = c(0, 25000)) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),
    legend.position = "right",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8)
  )
ggsave(filename = paste0(output_directory, name, "_Histogram_Cumulative_GeneExpression.pdf"), plot = cum_count, width = 22, height = 8, dpi = 300)

cum_count_plotly <- ggplotly(cum_count)
cum_count_plotly
saveWidget(cum_count_plotly, file = paste0(output_directory, name, "_Histogram_Cumulative_GeneExpression.html"))

# Calculate mean cumulative counts by condition
df_plot_grouped <- df_plot %>%
  group_by(bin_numeric, conds) %>%
  summarize(mean_cumulative_count = mean(cumulative_count))

# Plot
grouped_cum_count <- ggplot(df_plot_grouped, aes(x = bin_numeric, y = mean_cumulative_count, color = conds)) +
  geom_line() +
  scale_color_manual(values = color_dev[1:length(unique(df_plot_grouped$conds))]) +
  theme_minimal() +
  labs(
    title = "Histogram of Gene Expression Counts by Condition",
    x = "Bins (Size = 0.5)",
    y = "Mean Cumulative Count",
    color = "Condition"
  ) +
  coord_cartesian(ylim = c(0, 25000)) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),
    legend.position = "right",
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 10)
  )
ggsave(filename = paste0(output_directory, name, "_Grouped_Histogram_Cumulative_GeneExpression.pdf"), plot = grouped_cum_count, width = 16, height = 8, dpi = 300)

grouped_cum_count_plotly <- ggplotly(grouped_cum_count)
grouped_cum_count_plotly
saveWidget(grouped_cum_count_plotly, file = paste0(output_directory, name, "_Grouped_Histogram_Cumulative_GeneExpression.html"))

# Convert bins_matrix to a data frame
bins_df <- as.data.frame(bins_matrix)
bins_df$bin <- rownames(bins_matrix)

print(bins_df)

bins_long <- bins_df %>%
  pivot_longer(cols = -bin, names_to = "sample", values_to = "count")

# Extract the lower edge of each bin from the bin labels
# Ensure 'bin' is treated as a factor and ordered correctly
bins_long$bin <- factor(bins_long$bin, levels = unique(bins_long$bin), ordered = TRUE)

# Merge with conditions if available
bins_long <- merge(bins_long, conditions_ordered[, c("sample", "conds")], by = "sample")


# Plot using geom_col()
bin_count <- ggplot(bins_long, aes(x = bin, y = count, group = sample, color = sample)) +
  geom_line(aes(text = paste("Sample:", sample, "<br>Group:", conds, "<br>Bin:",bin,"<br>Count:",count))) +
  theme_minimal() +
  scale_color_manual(values = palette_colors) +
  labs(
    title = "Binned Histogram of Gene Expression Counts",
    x = "Bins (Size = 0.5)",
    y = "Gene Count",
    fill = "Sample"
  ) +
  coord_cartesian(ylim = c(0, 4000)) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    legend.position = "right",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8)
  )
ggsave(filename = paste0(output_directory, name, "_Binned_Histogram_GeneExpression.pdf"), plot = bin_count, width = 22, height = 8, dpi = 300)

bin_count_plotly <- ggplotly(bin_count, tooltip = "text")
bin_count_plotly
saveWidget(bin_count_plotly, file = paste0(output_directory, name, "_Binned_Histogram_GeneExpression.html"))

# Calculate mean counts per bin per condition
mean_bins <- bins_long %>%
  group_by(bin, conds) %>%
  summarize(mean_count = mean(count), .groups = "drop")

# Plot the mean histogram per condition
group_bin_count <- ggplot(mean_bins, aes(x = bin, y = mean_count, group = conds, color = conds)) +
  geom_line(aes(text = paste("Group:", conds, "<br>Bin:",bin,"<br>Averaged Count:",mean_count))) +
  theme_minimal() +
  scale_color_manual(values = palette_colors) +
  labs(
    title = "Mean Histogram of Gene Expression Counts per Condition",
    x = "Bins (Size = 0.5)",
    y = "Mean Gene Count",
    fill = "Condition"
  ) +
  coord_cartesian(ylim = c(0, 4000)) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    legend.position = "right",
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 10)
  )
ggsave(filename = paste0(output_directory, name, "_Grouped_Binned_Histogram_GeneExpression.pdf"), plot = group_bin_count, width = 16, height = 8, dpi = 300)

group_bin_count_plotly <- ggplotly(group_bin_count, tooltip = "text")
group_bin_count_plotly
saveWidget(group_bin_count_plotly, file = paste0(output_directory, name, "_Grouped_Binned_Histogram_GeneExpression.html"))











