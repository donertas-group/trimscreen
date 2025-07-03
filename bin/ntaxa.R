library(tidyverse)
library(scales)

wdir <- "/scratch/shire/ssd/pipeline/16s_nf_pipeline"

file.path(wdir, filename)

tab_mock12 = read.csv(file.path(wdir, "mock12/output/compare_runs/filtered_table.csv"))
tab_mock16 = read.csv(file.path(wdir, "mock16/output/compare_runs/filtered_table.csv"))
tab_mock21 = read.csv(file.path(wdir, "mock21/output/compare_runs/filtered_table.csv"))
tab_mock23 = read.csv(file.path(wdir, "mock23/output/compare_runs/filtered_table.csv"))

mock <- "Mock 12"; tab=tab_mock12
x <- "Phylum"; std.x=2
x <- "Family"; std.x=7
y <- "Genus"; std.y =11

mock <- "Mock 16"; tab=tab_mock16
x <- "Phylum"; std.x=20
x <- "Family"; std.x=43
y <- "Genus"; std.y =46

mock <- "Mock 21"; tab=tab_mock21
x <- "Phylum"; std.x=5
x <- "Family"; std.x=17
y <- "Genus"; std.y =17

mock <- "Mock 23"; tab=tab_mock23
x <- "Phylum"; std.x=5
x <- "Family"; std.x=17
y <- "Genus"; std.y =17

ggplot(tab, aes_string(x = x, y = y)) +
  geom_jitter(width = 0.1, height = 0.1, alpha = 0.4) +
  geom_vline(xintercept = std.x, color = "red", linetype = "dashed", size = 1) +
  geom_hline(yintercept = std.y, color = "red", linetype = "dashed", size = 1) +
  theme_minimal() +
  labs(title = paste(mock, ",", x, "vs", y),
       x = x, y = y) +
  #scale_x_continuous(labels = label_number(accuracy = 1)) +
  scale_y_continuous(labels = label_number(accuracy = 1)) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    panel.grid.minor = element_blank()
  )

# hist
ggplot(tab, aes_string(x = x)) +
  geom_histogram(binwidth = 1, fill = "skyblue", color = "black") +
  geom_vline(xintercept = std.x, color = "red", linetype = "dashed", size = 1) +
  labs(title = paste(mock, ",", x) ) +
  theme_minimal()


#######################################################
ds_mock16 =read.csv("../detection_summary.mock16.47family.txt",sep='\t')
ds_mock16 =read.csv("../detection_summary.mock16.51family.txt",sep='\t')
ds_mock16 =read.csv("../detection_summary.mock16.52family.txt",sep='\t')
ds_mock23 =read.csv("../detection_summary.mock23.20genus.txt",sep='\t')

mock="Mock-16"; ds = ds_mock16
mock="Mock-23"; ds = ds_mock23


ds_long <- ds %>%
  pivot_longer(cols = c(True_Detections, False_Detections, Missed_Detections), 
               names_to = "category", 
               values_to = "value")
ggplot(ds_long, aes(x = category, y = value, fill = category)) +
  geom_boxplot() +
  labs(title = paste(mock,", runs with mode richness, family"),
       x = "Category",
       y = "Value") +
  theme_minimal()
































