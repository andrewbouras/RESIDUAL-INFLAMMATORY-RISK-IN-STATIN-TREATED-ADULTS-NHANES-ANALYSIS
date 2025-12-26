# Install required R packages for RIR Study

packages <- c(
  "tidyverse",
  "survey",
  "tableone",
  "broom",
  "ggplot2",
  "here"
)

for (pkg in packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
    library(pkg, character.only = TRUE)
  }
}

cat("All required R packages installed.\n")

