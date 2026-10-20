#!/usr/bin/env Rscript
# 01_ingest_data.R  (v2)
# Phase 1: paired-cohort ingestion with local caching.
#
# Unit of ingestion = one platform (series matrix) OR one supplementary file.
# Multi-site series such as the SEQC releases yield many units, each mapped to
# HGNC symbols independently. Failures are per unit, never silent, and never
# replaced by synthetic data.
#
# Outputs in ./data/processed:
#   <tag>_expr.rds / <tag>_pheno.rds        per-unit gene x sample matrix + metadata
#   <tag>_expr_harmonized.rds               unit restricted to its group's gene index
#   harmonized_gene_index_<group>.rds       label-free gene index per cohort group
# and ./data/ingestion_manifest.csv (audit log).
# Unmappable files get a head-dump in ./data/diagnostics for manual inspection.

suppressPackageStartupMessages({
  library(GEOquery)
  library(Biobase)
})

DATA_DIR <- "./data"
RAW_DIR  <- file.path(DATA_DIR, "raw")
PROC_DIR <- file.path(DATA_DIR, "processed")
DIAG_DIR <- file.path(DATA_DIR, "diagnostics")
for (d in c(RAW_DIR, PROC_DIR, DIAG_DIR)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
options(timeout = 3600)

ACCESSIONS <- c("GSE56457", "GSE47774", "GSE135134", "GSE70353", "GSE211567")

# Cohort groups: each group gets its own gene index (label-free).
GROUPS <- list(
  development = c("GSE56457", "GSE47774"),
  validation  = c("GSE135134", "GSE70353"),
  stress      = c("GSE211567")
)

# Units excluded from the intersection that defines a group's gene index
# (e.g. targeted qPCR panels with ~1000 genes would collapse the gene space).
EXCLUDE_FROM_INDEX <- "qPCR"

# ---- generic helpers -------------------------------------------------------

strip_ansi <- function(x) gsub("\033\\[[0-9;]*m", "", x)

with_retry <- function(fn, tries = 4, wait = 10) {
  last <- NULL
  for (i in seq_len(tries)) {
    res <- tryCatch(fn(), error = function(e) e)
    if (!inherits(res, "error")) return(res)
    last <- res
    message("  attempt ", i, "/", tries, " failed: ", strip_ansi(conditionMessage(res)))
    if (i < tries) Sys.sleep(wait * i)
  }
  stop(strip_ansi(conditionMessage(last)), call. = FALSE)
}

sanitize <- function(x) gsub("[^A-Za-z0-9]+", "_", x)

find_symbol_col <- function(fd) {
  cands <- c("Gene Symbol", "Gene symbol", "GENE_SYMBOL", "Symbol", "SYMBOL",
             "gene_assignment", "Gene_Symbol", "ILMN_Gene", "GeneSymbol")
  hit <- intersect(cands, colnames(fd))
  if (length(hit)) return(hit[1])
  hit <- grep("symbol", colnames(fd), ignore.case = TRUE, value = TRUE)
  if (length(hit)) return(hit[1])
  NA_character_
}

clean_symbols <- function(x) {
  x <- as.character(x)
  x <- sub("\\s*(///|//|;|,).*$", "", x)
  x <- trimws(x)
  x[x %in% c("", "---", "NA", "NULL", "null")] <- NA_character_
  x
}


                    # Affymetrix ST arrays store "NM_xxx // SYMBOL // description /// ...".
extract_symbols <- function(tb, col) {
  x <- as.character(tb[[col]])
  if (col == "gene_assignment") {
    return(vapply(strsplit(x, " /// ", fixed = TRUE), function(s) {
      f <- strsplit(s[1], " // ", fixed = TRUE)[[1]]
      if (length(f) >= 2) trimws(f[2]) else NA_character_
    }, character(1)) |> (\(v) { v[v %in% c("", "---")] <- NA; v })())
  }
  clean_symbols(x)
}

collapse_probes <- function(mat, symbols) {
  keep <- !is.na(symbols) & symbols != ""
  mat <- mat[keep, , drop = FALSE]
  symbols <- symbols[keep]
  if (nrow(mat) == 0) stop("no rows left after symbol mapping")
  mu <- rowMeans(mat, na.rm = TRUE)
  ord <- order(symbols, -mu)
  mat <- mat[ord, , drop = FALSE]
  symbols <- symbols[ord]
  dup <- duplicated(symbols)
  mat <- mat[!dup, , drop = FALSE]
  rownames(mat) <- symbols[!dup]
  mat
}

need_pkgs <- function() {
  if (!requireNamespace("org.Hs.eg.db", quietly = TRUE) ||
      !requireNamespace("AnnotationDbi", quietly = TRUE))
    stop("Install Bioconductor packages 'org.Hs.eg.db' and 'AnnotationDbi'.")
}

org_map <- function(ids0, keytype) {
  need_pkgs()
  u <- unique(ids0)
  m <- suppressMessages(AnnotationDbi::mapIds(
    org.Hs.eg.db::org.Hs.eg.db, keys = u, column = "SYMBOL",
    keytype = keytype, multiVals = "first"))
  unname(m[ids0])
}


                    # ---- platform annotation (GPL) cache ---------------------------------------

.gpl_cache <- new.env()
get_gpl_table <- function(plat) {
  if (!is.null(.gpl_cache[[plat]])) return(.gpl_cache[[plat]])
  res <- tryCatch({
    gpl <- with_retry(function() getGEO(plat, destdir = RAW_DIR))
    tb <- Table(gpl)
    sc <- find_symbol_col(tb)
    if (is.na(sc) || !"ID" %in% colnames(tb)) NULL
    else data.frame(ID = as.character(tb$ID), symbol = extract_symbols(tb, sc),
                    stringsAsFactors = FALSE)
  }, error = function(e) NULL)
  .gpl_cache[[plat]] <- if (is.null(res)) FALSE else res
  .gpl_cache[[plat]]
}


                      # ---- supplementary file reading + ID mapping -------------------------------

read_supp_table <- function(f) {
  df <- read.delim(f, check.names = FALSE, stringsAsFactors = FALSE)
  auto_rn <- identical(rownames(df), as.character(seq_len(nrow(df))))
  if (auto_rn) {
    ids <- as.character(df[[1]]); df <- df[, -1, drop = FALSE]
  } else {
    ids <- rownames(df)
  }
  num <- vapply(df, is.numeric, logical(1))
  if (sum(num) < 2) stop("fewer than 2 numeric columns")
  list(ids = ids, mat = as.matrix(df[, num, drop = FALSE]),
       annot = df[, !num, drop = FALSE])
}


                      # Returns a symbol vector aligned to `ids`, plus how it was mapped.
map_supp_ids <- function(ids, annot, plats) {
  ids0 <- sub("\\.\\d+$", "", ids)
  frac <- function(p) mean(grepl(p, ids0))
  
  # 1. embedded symbol column
  sc <- find_symbol_col(annot)
  if (!is.na(sc)) return(list(sym = clean_symbols(annot[[sc]]), how = paste("column", sc)))
  
  # 2. Ensembl
  if (frac("^ENSG") > 0.5) return(list(sym = org_map(ids0, "ENSEMBL"), how = "Ensembl->org.Hs.eg.db"))
  
  # 3. probe IDs of any platform in this series (checked before Entrez, since
  #    some Affymetrix ST probeset IDs are purely numeric)
  best <- list(ov = 0)
  for (p in plats) {
    gt <- get_gpl_table(p)
    if (isFALSE(gt) || is.null(gt)) next
    ov <- mean(ids %in% gt$ID)
    if (ov > best$ov) best <- list(ov = ov, plat = p, gt = gt)
  }
  if (best$ov >= 0.5)
    return(list(sym = best$gt$symbol[match(ids, best$gt$ID)], how = paste("GPL", best$plat)))
 
  
  # 4. Entrez
  if (frac("^[0-9]+$") > 0.9) return(list(sym = org_map(ids0, "ENTREZID"), how = "Entrez->org.Hs.eg.db"))
 
  # 5. already symbols?
  need_pkgs()
  keys <- AnnotationDbi::keys(org.Hs.eg.db::org.Hs.eg.db, "SYMBOL")
  if (mean(ids0 %in% keys) > 0.5) return(list(sym = ids0, how = "symbols"))
  NULL
}


                      # Attach GEO metadata only when sample names match GSM ids or titles exactly.
attach_pheno <- function(sample_names, pheno_all, src) {
  base <- data.frame(sample_id = sample_names, source_file = src,
                     row.names = sample_names, stringsAsFactors = FALSE)
  if (is.null(pheno_all) || nrow(pheno_all) == 0) return(base)
  key <- ifelse(sample_names %in% rownames(pheno_all), sample_names,
                rownames(pheno_all)[match(sample_names, pheno_all$title)])
  if (anyNA(key)) return(base)
  out <- cbind(base, pheno_all[key, , drop = FALSE])
  rownames(out) <- sample_names
  out
}


                      # ---- unit builders ---------------------------------------------------------

unit_from_series <- function(es) {
  mat <- exprs(es); fd <- fData(es); plat <- annotation(es)
  sc <- find_symbol_col(fd)
  if (is.na(sc)) {
    gt <- get_gpl_table(plat)
    if (isFALSE(gt) || is.null(gt)) stop("no gene-symbol annotation for platform ", plat)
    sym <- gt$symbol[match(rownames(mat), gt$ID)]
    how <- paste("GPL", plat)
  } else {
    sym <- extract_symbols(fd, sc); how <- paste("fData column", sc)
  }
  list(expr = collapse_probes(mat, sym), pheno = pData(es),
       platform = plat, source = paste("series matrix;", how))
}

units_from_supp <- function(gse, plats, pheno_all) {
  files <- with_retry(function() getGEOSuppFiles(gse, baseDir = RAW_DIR, makeDirectory = TRUE))
  paths <- rownames(files)
  paths <- paths[grepl("\\.(txt|tsv|csv)(\\.gz)?$", paths, ignore.case = TRUE)]
  message("  ", length(paths), " tabular supplementary file(s)")
  out <- list(); fails <- list()
  for (f in paths) {
    nm <- basename(f)
    tag <- sanitize(sub("\\.(txt|tsv|csv)(\\.gz)?$", "", sub(paste0("^", gse, "_?"), "", nm)))
    res <- tryCatch({
      tb <- read_supp_table(f)
      mp <- map_supp_ids(tb$ids, tb$annot, plats)
      if (is.null(mp)) {
        writeLines(c(paste("File:", nm),
                     capture.output(print(utils::head(data.frame(id = tb$ids, tb$mat[, 1:min(4, ncol(tb$mat)), drop = FALSE]), 8)))),
                   file.path(DIAG_DIR, paste0(gse, "_", tag, "_head.txt")))
        stop("could not identify ID type (see data/diagnostics)")
      }
      m <- tb$mat; rownames(m) <- NULL
      m <- collapse_probes(m, mp$sym)
      list(expr = m, pheno = attach_pheno(colnames(m), pheno_all, nm),
           platform = NA_character_, source = paste(nm, ";", mp$how))
    }, error = function(e) e)
    if (inherits(res, "error")) fails[[tag]] <- strip_ansi(conditionMessage(res))
    else out[[tag]] <- res
  }
  list(units = out, fails = fails)
}


                    
                    # ---- per-accession ingestion ----------------------------------------------

ingest_one <- function(gse) {
  status_file <- file.path(PROC_DIR, paste0(gse, "_status.rds"))
  if (file.exists(status_file)) {
    message("[", gse, "] cached, skipping")
    return(readRDS(status_file))
  }
  message("[", gse, "] downloading / loading from cache")
  gl <- with_retry(function() getGEO(gse, GSEMatrix = TRUE, getGPL = FALSE, destdir = RAW_DIR))
  plats <- vapply(gl, annotation, character(1))
  pheno_all <- tryCatch(pData(gl[[1]]), error = function(e) NULL)
  empty <- vapply(gl, function(es) nrow(exprs(es)) == 0, logical(1))

  units <- list(); fails <- list()
  for (i in seq_along(gl)) {
    if (empty[i]) next
    tag <- plats[i]
    r <- tryCatch(unit_from_series(gl[[i]]), error = function(e) e)
    if (inherits(r, "error")) fails[[tag]] <- strip_ansi(conditionMessage(r))
    else units[[tag]] <- r
  }
  if (any(empty)) {
    message("  ", sum(empty), " empty series matrix(es); reading supplementary files")
    sp <- tryCatch(units_from_supp(gse, plats, pheno_all), error = function(e) e)
    if (inherits(sp, "error")) fails[["supplementary"]] <- strip_ansi(conditionMessage(sp))
    else { units <- c(units, sp$units); fails <- c(fails, sp$fails) }
  }

  multi <- length(units) + length(fails) > 1
  rows <- list()
  for (tag in names(units)) {
    u <- units[[tag]]
    full <- if (multi) paste0(gse, "__", sanitize(tag)) else gse
    saveRDS(u$expr,  file.path(PROC_DIR, paste0(full, "_expr.rds")))
    saveRDS(u$pheno, file.path(PROC_DIR, paste0(full, "_pheno.rds")))
    rows[[full]] <- data.frame(accession = full, parent = gse, platform = u$platform,
                               source = u$source, n_genes = nrow(u$expr),
                               n_samples = ncol(u$expr), status = "ok", note = "",
                               stringsAsFactors = FALSE)
  }
  for (tag in names(fails)) {
    full <- paste0(gse, "__", sanitize(tag))
    rows[[full]] <- data.frame(accession = full, parent = gse, platform = NA,
                               source = NA, n_genes = NA, n_samples = NA,
                               status = "failed", note = fails[[tag]],
                               stringsAsFactors = FALSE)
  }
  if (!length(rows)) stop("no units produced for ", gse)
  status <- do.call(rbind, rows); rownames(status) <- NULL
  
                   
                   # Cache only if everything succeeded, so failures are retried next run.
  if (all(status$status == "ok")) saveRDS(status, status_file)
  status
}

safe_ingest <- function(gse) {
  tryCatch(ingest_one(gse), error = function(e) {
    msg <- strip_ansi(conditionMessage(e))
    message("[", gse, "] FAILED: ", msg)
    data.frame(accession = gse, parent = gse, platform = NA, source = NA,
               n_genes = NA, n_samples = NA, status = "failed", note = msg,
               stringsAsFactors = FALSE)
  })
}


                   # ---- harmonization (per group, label-free) ---------------------------------

harmonize_groups <- function(manifest) {
  idx <- list()
  for (g in names(GROUPS)) {
    parents <- GROUPS[[g]]
    ok <- manifest[manifest$status == "ok" & manifest$parent %in% parents, ]
    missing <- setdiff(parents, ok$parent)
    if (length(missing)) {
      warning("Group '", g, "' skipped: no usable data for ", paste(missing, collapse = ", "))
      next
    }
    use <- ok$accession[!grepl(EXCLUDE_FROM_INDEX, ok$accession)]
    gene_sets <- lapply(use, function(a)
      rownames(readRDS(file.path(PROC_DIR, paste0(a, "_expr.rds")))))
    shared <- sort(Reduce(intersect, gene_sets))
    message("Group '", g, "': ", length(shared), " genes shared across ", length(use), " unit(s)")
    if (!length(shared)) { warning("Group '", g, "' has an empty gene intersection"); next }
    saveRDS(shared, file.path(PROC_DIR, paste0("harmonized_gene_index_", g, ".rds")))
    for (a in ok$accession) {
      m <- readRDS(file.path(PROC_DIR, paste0(a, "_expr.rds")))
      sub <- m[intersect(shared, rownames(m)), , drop = FALSE]
      if (nrow(sub) < length(shared))
        message("  note: ", a, " covers only ", nrow(sub), "/", length(shared), " index genes")
      saveRDS(sub, file.path(PROC_DIR, paste0(a, "_expr_harmonized.rds")))
    }
    idx[[g]] <- shared
  }
  invisible(idx)
}


                        # ---- sanity checks ---------------------------------------------------------

run_checks <- function(manifest) {
  m <- matrix(c(1, 2, 10, 20, 5, 6), nrow = 3, byrow = TRUE,
              dimnames = list(c("p1", "p2", "p3"), c("s1", "s2")))
  out <- collapse_probes(m, c("A", "A", NA))
  stopifnot(identical(rownames(out), "A"), all(out["A", ] == c(10, 20)))
  stopifnot(identical(clean_symbols(c("TP53 /// X", "---", " BRCA1 ")), c("TP53", NA, "BRCA1")))
  tb <- data.frame(gene_assignment = c("NM_1 // GENE1 // desc /// NR_2 // GENE9 // d", "---"))
  stopifnot(identical(extract_symbols(tb, "gene_assignment"), c("GENE1", NA)))
  for (a in manifest$accession[manifest$status == "ok"]) {
    e <- readRDS(file.path(PROC_DIR, paste0(a, "_expr.rds")))
    p <- readRDS(file.path(PROC_DIR, paste0(a, "_pheno.rds")))
    stopifnot(!anyDuplicated(rownames(e)), !anyNA(rownames(e)))
    if (ncol(e) != nrow(p))
      warning(a, ": expression columns (", ncol(e), ") != metadata rows (", nrow(p), ")")
  }
  message("All sanity checks passed")
}


                        # ---- main ------------------------------------------------------------------

if (sys.nframe() == 0) {
  manifest <- do.call(rbind, lapply(ACCESSIONS, safe_ingest))
  rownames(manifest) <- NULL
  write.csv(manifest, file.path(DATA_DIR, "ingestion_manifest.csv"), row.names = FALSE)
  print(manifest[, c("accession", "platform", "n_genes", "n_samples", "status")])
  if (any(manifest$status == "failed")) {
    warning("Some units failed; read ingestion_manifest.csv and data/diagnostics before Phase 2.")
  }
  if (any(manifest$status == "ok")) {
    harmonize_groups(manifest)
    run_checks(manifest)
  }
}
