# terminal command: `Rscript main_script.R`

# check if packages are installed and install if not

packages <- c("tidyverse", "knitr", "kableExtra", "nlme",
              "emmeans", "standardize", "gridExtra", "effectsize",
              "nlstools", "htmltools", "rmarkdown")

for(pckg in packages) {
  if(!pckg %in% rownames(installed.packages())) {
    install.packages(pckg)
  }
}


# render RMarkdown files

rmarkdown::render(
  input = "analysis/lmm/LMM_AM_byLayer_AllTrials_KIT.Rmd",  
  output_file = "analysis/lmm/LMM_AM_byLayer_AllTrials_KIT.html")

rmarkdown::render(
  input = "analysis/lmm/LMM_CL_byLayer_AllTrials_KIT.Rmd",  
  output_file = "analysis/lmm/LMM_CL_byLayer_AllTrials_KIT")