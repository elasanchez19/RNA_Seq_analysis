


load_config <- function(path = "config/config.yaml") yaml::read_yaml(path)
ensure_dir <- funtion(...) dir.create(file.path(...), recursive = TRUE,showWarnings = FALSE)

pal_condition <- c(Sham = "#7f8c8d", SNI = "#d35400")

theme_repo <- function(base_size = 12) {
  theme_bw(base_size = base_size) + 
    theme(panel.grid.minor = element_blank(),
          strip.background = element_rect(fill = "grey92", colour = NA))
}

# ----- Counts lecture -----
# admit metadata with column "counts_file" + folder data.counts_path with
# a file per sample (col 1 = gen, col 2 = counts; .csv/.csv.gz, with ot without header)

detect_sep <- function(f) {
  first <- readLines(f, n=1, warn = FALSE)
  if (length(first) == 0) stop("Empty file: ", basename(f))
  if (grepl("\t", first, fixed = TRUE)) "\t"
  else if (grepl(",", first, fixed = TRUE)) ","
  else stop("Cannot detect separator in: " , basename(f))
}

read_count_file <- function(f){
  sep <- detect_sep(f)
  d <- read.table(f, header=FALSE, sep=sep, stringsAsFactors = FALSE,
                  quote = "\"", comment.char = "", check.names = FALSE, fill=TRUE)
  if (ncol(d) < 2) stop("File", basename(f), "has fewer than columns")
  if (ncol(d) > 2)
    warning(basename(f), " has ", ncol(d), " columns; using columns 1 and 2 as gene and count.")
  
  if (is.na(suppressWarnings(as.numeric(d[1,2])))) d <- d[-1, drop = FALSE]
  
  counts <- suppressWarnings(as.numeric(d[[2]]))
  if (anyNA(counts)) stop("Non-numeric counts found in: ", basename(f))
  
  out <- data.frame(gene = d[[1]], count = counts, stringsAsFactors = FALSE)
  
  #remove   featurecounts summary categories, not genes
  out <- out[!grepl("^__", out$gene), ]
  
  if (anyDuplicated(out$gene)) {
    warning(basename(f), ": duplicated gene names, counts were summed.")
    agg <- rowsum(out$count, out$gene)
    out <- data.frame(gene = rownames(agg), count = as.numeric(agg[,1]),
                      stringsAsFactors = FALSE)
  }
  out
}

read_counts_from_files <- function(dir, files, ids) {
  paths <- file.path(dir, files)
  bad <- !file.exists(paths)
  if (any(bad)) stop("There are not these files in ", dir, ":\n  ",
                     paste(files[bad], collapse= "\n  "),
                     "\nFolder content:\n  ", paste(head(list.files(dir), 20), collapse = "\n  "))
  lst <- Map(function(f, id) { x <- read_count_file(f); names(x)[2] <- id; x }, paths, ids)
  m <- Reduce(function(a, b) merge(a, b, by = "gene", all = TRUE), lst)
  stopifnot(!anyDuplicated(m$gene))
  rownames(m) <- m$gene; m$gene <- NULL
  m[is.na(m)] <-0
  as.matrix(m)
}


# ----- Genes annotation -----
get_annotation <- function(gene_ids, cfg) {
  if (cfg$data$gene_id_type == "symbol") {
    return(data.frame(gene_id = gene_ids, symbol = gene_ids, stringsAsFactors = FALSE))
  }
  suppressPackageStartupMessages(library(cfg$data$organism_db, character.only = TRUE))
  db <- get(cfg$data$organism_db)
  clean <- sub("\\..*$", "", gene_ids)
  sym <- AnnotationDbi::mapIds(db, keys = clean, keytype = "ENSEMBL",
                               column = "SYMBOL", multiVals = "first")
  data.frame(gene_id = gene_ids, symbol = ifelse(is.na(sym), gene_ids, sym),
             stringsAsFactors = FALSE)
}

# ----- Metadata -----

normalize_meta <- function(meta) {
  sx <- toupper(substr(trimws(meta$sex),1,1))
  meta$sex <- sx
  cn <- toupper(trimws(meta$condition))
  meta$condition <- ifelse(cn == "SHAM", "Sham", ifelse(cn == "SNI", "SNI", meta$condition))
  meta
}

load_data <- function(cfg) {
 #counts tells you what was measured.
 #meta tells you under which conditions each measurement was taken.
 #annot tells you which biological gene each row represents.
  meta <- normalize_meta(read.csv(cfg$data$metadata_path, stringsAsFactors = F))
  stopifnot(all(c("sample_id", "sex", "condition") %in% names(meta)))
  if ("counts_file" %in% names(meta) && dir.exists(cfg$data$counts_path)) {
    # for separate samples_counts.cvs files
    counts <- read_counts_from_files(cfg$data$counts_path, meta$counts_file, meta$sample_id)
  } else {
    stop("There is not the column count_file", cf$data$counts_path)
  }
  meta$condition <- factor(meta$condition, levels=cfg$factors$condition_levels)
  meta$sex <- factor(meta$sex, levels = cfg$factors$sex_levels)
  meta$group <- interaction(meta$sex, meta$condition, sep="_")
  rownames(meta) <- meta$sample_id
  list(counts = counts, meta=meta, annot=get_annotation(rownames(counts), cfg))
}


