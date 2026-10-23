source("scripts/utils.R")

install.packages("BiocManager")
BiocManager::install("edgeR")
library(edgeR)
library(ggplot2)

cfg <- load_config()
d <- load_data(cfg)

print(table(sex = d$meta$sex, condition = d$meta$condition))

# filter genes

counts <- d$counts; meta <- d$meta; annot <- d$annot
ensure_dir(cfg$paths$processed)


# Size libraries

lib <- data.frame(sample= colnames(counts),
                  lib_size_M = colSums(counts) / 1e6,
                  sex = meta[colnames(counts), "sex"],
                  condition = meta[colnames(counts), "condition"])
p_lib <- ggplot(lib, aes(reorder(sample, lib_size_M), lib_size_M, fill = condition)) +
  geom_col() + coord_flip() + facet_wrap(~ sex, scales = "free_y") +
  scale_fill_manual(values = pal_condition) + 
  labs(x = NULL, y = "Library size (millions of counts)", title = "Sequencing depth per sample") + 
  theme_repo()
p_lib
save_plot(p_lib, "library_sizes", cfg, 8, 4)
write_table(lib, "library_size.csv", cfg)



