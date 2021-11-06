# terminal command: `Rscript package_install.R`
# or run the code directly 
packages <- c("tidyverse", "ggpubr", "rstatix", "knitr", "kableExtra", "nlme",
              "emmeans", "standardize", "sjPlot", "gridExtra", "effectsize",
              "plotly", "nlstools", "htmltools")

for(pack in packages) {
  if(!pack %in% rownames(installed.packages())) {
    install.packages(pack)
  }
}
