pkgs <- c("tidyverse", "cowplot", "patchwork", "viridis", "RColorBrewer",
          "wesanderson", "ggdist", "ggridges", "scales", "ggExtra",
          "ggpubr", "tidyplots")

to_install <- pkgs[!pkgs %in% installed.packages()[, "Package"]]
if (length(to_install) > 0) {
  cat("Installing:", paste(to_install, collapse = ", "), "\n")
  install.packages(to_install, repos = "https://cran.r-project.org")
} else {
  cat("All packages already installed.\n")
}

# ComplexHeatmap via Bioconductor
if (!"ComplexHeatmap" %in% installed.packages()[, "Package"]) {
  BiocManager::install("ComplexHeatmap", update = FALSE, ask = FALSE)
}

cat("\nDone.\n")
