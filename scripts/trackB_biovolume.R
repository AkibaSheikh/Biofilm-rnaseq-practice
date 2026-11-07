
library(tidyverse)

m  <- read.csv("results/trackB_measurements.csv")
dz <- 0.4    # slice spacing in micrometres, from the dataset's metadata.mat

m$depth_um <- (seq_len(nrow(m)) - 1) * dz
biovolume_um3 <- sum(m$Area) * dz
biovolume_um3   # should be about 45

p <- ggplot(m, aes(depth_um, Area)) + geom_line() +
  labs(x = "Depth (um)", y = "Thresholded area (um^2)",
       title = "Biofilm area per slice") +
  theme_bw()

ggsave("figures/trackB_area_vs_depth.png", p, width = 6, height = 4, dpi = 300)
