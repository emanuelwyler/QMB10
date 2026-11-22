
install.packages(c(
"Seurat",
"dplyr",
"ggplot2",
"tibble",
"tidyr",
"igraph",
"reshape2",
"dendextend",
"stringr",
"gridExtra",
"cowplot",
"data.table",
"pals",
"presto",
"purrr",
"ggpubr",
"RColorBrewer",
"lme4"
))

install.packages(
  "BiocManager",
  repos = "https://cloud.r-project.org"
)

BiocManager::install(
  c("DESeq2", "apeglm"),
  ask = FALSE,
  update = FALSE
)

install.packages('VoltRon', repos = c('https://bimsbbioinfo.r-universe.dev',
                                      'https://bioc.r-universe.dev',
                                      'https://cloud.r-project.org'))


if (!requireNamespace("VoltRon", quietly = TRUE)) {
  stop("VoltRon could not be installed or loaded!")
}