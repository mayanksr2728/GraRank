# 03_grassmann_engine.R -- Phase 3: Grassmann Manifold Subspace Alignment & Benchmarking (Fully Fixed)

need <- c("glmnet", "e1071", "pROC")
miss <- need[!vapply(need, requireNamespace, logical(1), quietly = TRUE)]
if (length(miss)) stop("Install first: install.packages(c(", paste0('"', miss, '"', collapse = ", "), "))")

PROC_DIR <- "./data/processed"
RES_DIR <- "./data/results"
dir.create(RES_DIR, recursive = TRUE, showWarnings = FALSE)
set.seed(2026)

---- 1. Grassmann Subspace Projection Engine --------------------------------

project_grassmann <- function(X_tr, X_te, k = 20) {
  svd_tr <- svd(t(X_tr), nu = min(k, nrow(X_tr) - 1, ncol(X_tr) - 1), nv = 0)
  U_k <- svd_tr$u
  Z_tr <- X_tr %% U_k
  Z_te <- X_te %% U_k
  colnames(Z_tr) <- paste0("PC", seq_len(ncol(Z_tr)))
  colnames(Z_te) <- paste0("PC", seq_len(ncol(Z_te)))
  list(Z_tr = Z_tr, Z_te = Z_te, U_k = U_k)
}

fit_svm_subspace <- function(Z_tr, y_tr) {
  m <- e1071::svm(Z_tr, factor(y_tr, levels = c(0, 1)), kernel = "linear", cost = 1, scale = TRUE)
  dv_tr <- as.numeric(attr(predict(m, Z_tr, decision.values = TRUE), "decision.values"))
  list(m = m, flip = mean(dv_tr[y_tr == 1]) < mean(dv_tr[y_tr == 0]))
}

pred_svm_subspace <- function(fit, Z_te) {
  p <- predict(fit$m, Z_te, decision.values = TRUE)
  dv <- as.numeric(attr(p, "decision.values"))
  if (fit$flip) -dv else dv
}

rank_transform <- function(m) {
  if (anyNA(m)) stop("NA values present; cannot rank")
  r <- apply(m, 2, function(x) rank(x, ties.method = "average"))
  if (!is.matrix(r)) {
    r <- matrix(r, nrow = nrow(m), ncol = ncol(m))
  }
  r <- r / nrow(m)
  dimnames(r) <- dimnames(m)
  r
}

top_genes_t <- function(X, y, n) {
  X1 <- X[y == 1, , drop = FALSE]
  X0 <- X[y == 0, , drop = FALSE]
  m1 <- colMeans(X1)
  m0 <- colMeans(X0)
  v1 <- (colSums(X1^2) - nrow(X1) * m1^2) / (nrow(X1) - 1)
  v0 <- (colSums(X0^2) - nrow(X0) * m0^2) / (nrow(X0) - 1)
  t <- (m1 - m0) / sqrt(pmax(v1, 0) / nrow(X1) + pmax(v0, 0) / nrow(X0) + 1e-12)
  names(sort(abs(t), decreasing = TRUE))[seq_len(min(n, ncol(X)))]
}

fit_lasso <- function(X, y) {
  nf <- min(5, min(table(y)))
  if (nf < 3) stop("too few samples per class for cross-validation")
  glmnet::cv.glmnet(X, factor(y, levels = c(0, 1)), family = "binomial", alpha = 1, nfolds = nf, type.measure = "deviance")
}

pred_lasso <- function(cv, X) {
  as.numeric(predict(cv, newx = X, s = "lambda.min", type = "link"))
}

safe_auc <- function(y, s) {
  r <- pROC::roc(y, s, levels = c(0, 1), direction = "<", quiet = TRUE)
  ci <- suppressWarnings(as.numeric(pROC::ci.auc(r, method = "delong")))
  list(auc = as.numeric(pROC::auc(r)), lo = ci[1], hi = ci[3])
}

run_grassmann_methods <- function(tr, te, k = 20) {
  proj <- project_grassmann(tr$X, te$X, k = k)
  res <- list()
  
  r_lasso <- tryCatch({
    fit_l <- fit_lasso(proj$Z_tr, tr$y)
    pred_l <- pred_lasso(fit_l, proj$Z_te)
    a_l <- safe_auc(te$y, pred_l)
    data.frame(method = "Grassmann+LASSO", auc = a_l$auc, auc_lo = a_l$lo, auc_hi = a_l$hi, accuracy = mean((pred_l > 0) == (te$y == 1)), note = "")
  }, error = function(e) {
    data.frame(method = "Grassmann+LASSO", auc = NA, auc_lo = NA, auc_hi = NA, accuracy = NA, note = conditionMessage(e))
  })
  
  r_svm <- tryCatch({
    fit_s <- fit_svm_subspace(proj$Z_tr, tr$y)
    pred_s <- pred_svm_subspace(fit_s, proj$Z_te)
    a_s <- safe_auc(te$y, pred_s)
    data.frame(method = "Grassmann+SVM", auc = a_s$auc, auc_lo = a_s$lo, auc_hi = a_s$hi, accuracy = mean((pred_s > 0) == (te$y == 1)), note = "")
  }, error = function(e) {
    data.frame(method = "Grassmann+SVM", auc = NA, auc_lo = NA, auc_hi = NA, accuracy = NA, note = conditionMessage(e))
  })
  
  rbind(r_lasso, r_svm)
}

manifest <- read.csv("./data/ingestion_manifest.csv", stringsAsFactors = FALSE)
ok <- manifest[manifest$status == "ok", ]

load_expr <- function(acc) readRDS(file.path(PROC_DIR, paste0(acc, "expr_harmonized.rds")))
load_index <- function(g) readRDS(file.path(PROC_DIR, paste0("harmonized_gene_index", g, ".rds")))
mk <- function(X, y, batch) list(X = X, y = as.integer(y), batch = rep(batch, length(y)))

seqc_letter <- function(x) {
  m <- regexpr("[A-D][0-9]+", x)
  ifelse(m > 0, substr(x, m + 1, m + 1), NA_character_)
}

build_seqc <- function(cap = 20) {
  genes <- load_index("development")
  take <- function(acc, cap) {
    m <- load_expr(acc)
    stopifnot(all(genes %in% rownames(m)))
    m <- m[genes, , drop = FALSE]
    lab <- seqc_letter(colnames(m))
    keep <- which(lab %in% c("A", "B"))
    if (!is.null(cap)) {
      keep <- unlist(lapply(c("A", "B"), function(l) {
        k <- keep[lab[keep] == l]
        if (length(k) > cap) sample(k, cap) else k
      }))
    }
    list(X = t(rank_transform(m[, keep, drop = FALSE])), y = as.integer(lab[keep] == "B"))
  }
  comb <- function(l) {
    list(X = do.call(rbind, lapply(l, function(z) z$X)), y = unlist(lapply(l, function(z) z$y)))
  }
  accs_arr <- ok$accession[ok$parent == "GSE56457" & !grepl("qPCR", ok$accession) & !grepl("GPL18522", ok$accession)]
  accs_rna <- ok$accession[ok$parent == "GSE47774"]
  arr <- comb(lapply(accs_arr, take, cap = NULL))
  rna <- comb(lapply(accs_rna, take, cap = cap))
  list(arr = mk(arr$X, arr$y, "array"), rna = mk(rna$X, rna$y, "rnaseq"))
}

build_metsim <- function() {
  genes <- load_index("validation")
  prep <- function(acc) {
    e <- load_expr(acc)
    ph <- readRDS(file.path(PROC_DIR, paste0(acc, "pheno.rds")))
    i <- match(colnames(e), rownames(ph))
    stopifnot(!anyNA(i))
    ph <- ph[i, , drop = FALSE]
    stopifnot(all(genes %in% rownames(e)))
    bmi <- as.numeric(ph[["bmi:ch1"]])
    y <- ifelse(bmi >= 30, 1L, ifelse(bmi < 25, 0L, NA_integer))
    keep <- which(!is.na(y))
    list(X = t(rank_transform(e[genes, keep, drop = FALSE])), y = y[keep], id = as.character(ph[["title"]][keep]))
  }
  res_rna <- prep("GSE135134")
  res_arr <- prep("GSE70353")
  list(rna = res_rna, arr = res_arr)
}

drop_overlap <- function(te, tr) {
  k <- !(te$id %in% tr$id)
  message(" removed ", sum(!k), " test individuals that are also in the training cohort")
  list(X = te$X[k, , drop = FALSE], y = te$y[k], id = te$id[k])
}

build_fever <- function() {
  ph <- readRDS(file.path(PROC_DIR, "GSE211567_pheno.rds"))
  e <- load_expr("GSE211567")
  genes <- load_index("stress")
  stopifnot(all(colnames(e) == rownames(ph)), all(genes %in% rownames(e)))
  lab <- ph[["infection type:ch1"]]
  site <- ph[["site:ch1"]]
  keep <- which(ph$geo_matched & lab %in% c("Bacterial", "Viral"))
  X <- t(rank_transform(e[genes, keep, drop = FALSE]))
  list(X = X, y = as.integer(lab[keep] == "Bacterial"), site = site[keep])
}

tasks <- list()
add_task <- function(name, direction, tr, te) {
  for (s in list(tr, te)) if (length(unique(s$y)) < 2) stop(name, " ", direction, ": a split has one class")
  tasks[[length(tasks) + 1]] <<- list(name = name, direction = direction, tr = tr, te = te)
}

sq <- build_seqc()
add_task("SEQC B-vs-A", "array -> RNA-seq", sq$arr, sq$rna)
add_task("SEQC B-vs-A", "RNA-seq -> array", sq$rna, sq$arr)

ms <- build_metsim()
te1 <- drop_overlap(ms$arr, ms$rna)
add_task("METSIM obese-vs-lean", "RNA-seq -> array", mk(ms$rna$X, ms$rna$y, "rnaseq"), mk(te1$X, te1$y, "array"))
te2 <- drop_overlap(ms$rna, ms$arr)
add_task("METSIM obese-vs-lean", "array -> RNA-seq", mk(ms$arr$X, ms$arr$y, "array"), mk(te2$X, te2$y, "rnaseq"))

fv <- build_fever()
for (s0 in c("Sri_Lanka", "United_States")) {
  a <- fv$site == s0
  add_task("Fever bacterial-vs-viral", paste(s0, "->", setdiff(unique(fv$site), s0)), mk(fv$X[a, , drop = FALSE], fv$y[a], s0), mk(fv$X[!a, , drop = FALSE], fv$y[!a], setdiff(unique(fv$site), s0)))
}

results <- do.call(rbind, lapply(tasks, function(t) {
  message("\n== Grassmann Engine | ", t$name, " | ", t$direction)
  cbind(task = t$name, direction = t$direction, n_train = length(t$tr$y), n_test = length(t$te$y), run_grassmann_methods(t$tr, t$te, k = 20))
}))

rownames(results) <- NULL
results[, c("auc", "auc_lo", "auc_hi", "accuracy")] <- round(results[, c("auc", "auc_lo", "auc_hi", "accuracy")], 3)

write.csv(results, file.path(RES_DIR, "phase3_grassmann.csv"), row.names = FALSE)
print(results, row.names = FALSE)
