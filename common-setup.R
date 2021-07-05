library(tidyverse, warn.conflicts = FALSE)
library(patchwork)
library(ggplot2)
library(patchwork)

here::i_am("common-setup.R")

if (!"ukcovidtools" %in% c(devtools::dev_packages(),rownames(installed.packages()))) {
  if (fs::dir_exists("~/Git/uk-covid-datatools/")) {
    devtools::load_all("~/Git/uk-covid-datatools/")
  } else {
    devtools::install_github("terminological/uk-covid-datatools")
  }
}
ukcovidtools::setup(here::here("config.yml"))

reset = function() {
  devtools::load_all("~/Git/uk-covid-datatools/")
  ukcovidtools::reload(here::here("config.yml"))
}

if (!"standardPrintOutput" %in% c(devtools::dev_packages(),rownames(installed.packages()))) {
  if (fs::dir_exists("~/Git/standard-print-output/")) {
    devtools::load_all("~/Git/standard-print-output/")
  } else {
    devtools::install_github("terminological/standard-print-output")
  }
}
standardPrintOutput::setDefaults()



