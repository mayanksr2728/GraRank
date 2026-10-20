install.packages(c("glmnet", "e1071", "pROC"))
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("sva")   # optional: enables the ComBat baseline


# 02_rank_transforms.R -- Phase 2: within-sample ranks + baseline classifiers
# Tasks (class 1 listed first):
#   SEQC     : B vs A      array <-> RNA-seq
#   METSIM   : obese (BMI>=30) vs lean (BMI<25), RNA-seq <-> array,
#              individuals present in the training cohort removed from the test cohort
#   Fever    : Bacterial vs Viral, Sri_Lanka <-> United_States
# Ranks are computed per sample, so no information crosses samples or folds.

need <- c("glmnet", "e1071", "pROC")
miss <- need[!vapply(need, requireNamespace, logical(1), quietly = TRUE)]
if (length(miss)) stop("Install first: install.packages(c(",
                       paste0('"', miss, '"', collapse = ", "), "))")
HAVE_SVA <- requireNamespace("sva", quietly = TRUE)
if (!HAVE_SVA) message("sva not installed: ComBat baseline will be skipped")

PROC_DIR <- "./data/processed"
RES_DIR  <- "./data/results"
dir.create(RES_DIR, recursive = TRUE, showWarnings = FALSE)
set.seed(2026)

# ---- core functions --------------------------------------------------------
# m: genes x samples. Returns relative ranks in (0, 1], ties averaged.
rank_transform <- function(m) {
  if (anyNA(m)) stop("NA values present; cannot rank")
  r <- apply(m, 2, rank, ties.method = "average")
  r <- r / nrow(m)
  dimnames(r) <- dimnames(m)
  r
}

top_genes_t <- function(X, y, n) {
  X1 <- X[y == 1, , drop = FALSE]; X0 <- X[y == 0, , drop = FALSE]
  m1 <- colMeans(X1); m0 <- colMeans(X0)
  v1 <- (colSums(X1^2) - nrow(X1) * m1^2) / (nrow(X1) - 1)
  v0 <- (colSums(X0^2) - nrow(X0) * m0^2) / (nrow(X0) - 1)
  t <- (m1 - m0) / sqrt(pmax(v1, 0) / nrow(X1) + pmax(v0, 0) / nrow(X0) + 1e-12)
  names(sort(abs(t), decreasing = TRUE))[seq_len(min(n, ncol(X)))]
}

# ---- classifiers (all fitted on training data only) ------------------------
fit_lasso <- function(X, y) {
  nf <- min(5, min(table(y)))
  if (nf < 3) stop("too few samples per class for cross-validation")
  glmnet::cv.glmnet(X, factor(y, levels = c(0, 1)), family = "binomial",
                    alpha = 1, nfolds = nf, type.measure = "deviance")
}
pred_lasso <- function(cv, X) as.numeric(predict(cv, newx = X, s = "lambda.min", type = "response"))

fit_svm <- function(X, y, ntop = 1000) {
  genes <- top_genes_t(X, y, ntop)
  m <- e1071::svm(X[, genes, drop = FALSE], factor(y, levels = c(0, 1)),
                  kernel = "linear", cost = 1, scale = TRUE)
  list(m = m, genes = genes)
}
pred_svm <- function(fit, X) {
  p <- predict(fit$m, X[, fit$genes, drop = FALSE], decision.values = TRUE)
  dv <- as.numeric(attr(p, "decision.values"))
  r <- suppressWarnings(cor(dv, as.numeric(as.character(p))))
  if (!is.na(r) && r < 0) dv <- -dv    # orient so that higher = class 1
  dv
}

ktsp_fit <- function(X, y, k = 5, ncand = 200) {
  cand <- top_genes_t(X, y, ncand)
  Xc <- X[, cand, drop = FALSE]
  X1 <- Xc[y == 1, , drop = FALSE]; X0 <- Xc[y == 0, , drop = FALSE]
  g <- ncol(Xc)
  P1 <- matrix(0, g, g); P0 <- matrix(0, g, g)
  for (i in seq_len(g)) {
    P1[i, ] <- colMeans(X1[, i] < X1)    # P(x_i < x_j | class 1)
    P0[i, ] <- colMeans(X0[, i] < X0)
  }
  D <- P1 - P0
  S <- abs(D); S[lower.tri(S, diag = TRUE)] <- 0
  pairs <- list(); used <- character(0)
  for (o in order(S, decreasing = TRUE)) {
    if (length(pairs) >= k || S[o] <= 0) break
    i <- (o - 1) %% g + 1; j <- (o - 1) %/% g + 1
    if (cand[i] %in% used || cand[j] %in% used) next
    pairs[[length(pairs) + 1]] <- list(i = cand[i], j = cand[j], s = sign(D[i, j]))
    used <- c(used, cand[i], cand[j])
  }
  list(pairs = pairs)
}
ktsp_predict <- function(fit, X) {
  if (!length(fit$pairs)) return(rep(0, nrow(X)))
  votes <- vapply(fit$pairs, function(p) {
    lt <- X[, p$i] < X[, p$j]
    ifelse(lt, p$s, -p$s)
  }, numeric(nrow(X)))
  votes <- matrix(votes, nrow = nrow(X))
  rowSums(votes) / length(fit$pairs)
}

# ComBat (batch = platform/site) on train+test ranks, then LASSO on corrected data.
# Uses test-set feature distributions (no labels) -> transductive; flagged in output.
combat_lasso <- function(tr, te) {
  dat <- t(rbind(tr$X, te$X))
  batch <- c(tr$batch, te$batch)
  keep <- rep(TRUE, nrow(dat))
  for (b in unique(batch)) {
    sub <- dat[, batch == b, drop = FALSE]
    v <- rowSums((sub - rowMeans(sub))^2) / (ncol(sub) - 1)
    keep <- keep & v > 1e-10
  }
  adj <- sva::ComBat(dat[keep, ], batch = batch, par.prior = TRUE, prior.plots = FALSE)
  ntr <- nrow(tr$X)
  Xtr <- t(adj[, seq_len(ntr)]); Xte <- t(adj[, ntr + seq_len(nrow(te$X))])
  pred_lasso(fit_lasso(Xtr, tr$y), Xte)
}

# ---- evaluation ------------------------------------------------------------
safe_auc <- function(y, s) {
  r <- pROC::roc(y, s, levels = c(0, 1), direction = "<", quiet = TRUE)
  ci <- suppressWarnings(as.numeric(pROC::ci.auc(r, method = "delong")))
  list(auc = as.numeric(pROC::auc(r)), lo = ci[1], hi = ci[3])
}

run_methods <- function(tr, te) {
  res <- list()
  add <- function(method, f, thr, note = "") {
    r <- tryCatch({
      s <- f(); a <- safe_auc(te$y, s)
      data.frame(method = method, auc = a$auc, auc_lo = a$lo, auc_hi = a$hi,
                 accuracy = mean((s > thr) == (te$y == 1)), note = note)
    }, error = function(e)
      data.frame(method = method, auc = NA, auc_lo = NA, auc_hi = NA,
                 accuracy = NA, note = conditionMessage(e)))
    res[[length(res) + 1]] <<- r
  }
  add("LASSO",  function() pred_lasso(fit_lasso(tr$X, tr$y), te$X), 0.5)
  add("SVM",    function() pred_svm(fit_svm(tr$X, tr$y), te$X), 0)
  add("k-TSP",  function() ktsp_predict(ktsp_fit(tr$X, tr$y), te$X), 0)
  if (HAVE_SVA)
    add("ComBat+LASSO", function() combat_lasso(tr, te), 0.5, "transductive (uses test features)")
  do.call(rbind, res)
}

# ---- data loading ----------------------------------------------------------
manifest <- read.csv("./data/ingestion_manifest.csv", stringsAsFactors = FALSE)
ok <- manifest[manifest$status == "ok", ]
load_expr <- function(acc) readRDS(file.path(PROC_DIR, paste0(acc, "_expr_harmonized.rds")))
load_index <- function(g) readRDS(file.path(PROC_DIR, paste0("harmonized_gene_index_", g, ".rds")))
mk <- function(X, y, batch) list(X = X, y = as.integer(y), batch = rep(batch, length(y)))

seqc_letter <- function(x) {
  m <- regexpr("_[A-D]_[0-9]+", x)
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
    if (!is.null(cap))
      keep <- unlist(lapply(c("A", "B"), function(l) {
        k <- keep[lab[keep] == l]
        if (length(k) > cap) sample(k, cap) else k
      }))
    message(sprintf("  %-40s A/B samples used: %d", acc, length(keep)))
    list(X = t(rank_transform(m[, keep, drop = FALSE])), y = as.integer(lab[keep] == "B"))
  }
  comb <- function(l) list(X = do.call(rbind, lapply(l, function(z) z$X)),
                           y = unlist(lapply(l, function(z) z$y)))
  message("SEQC units:")
  arr <- comb(lapply(ok$accession[ok$parent == "GSE56457" &
                                  !grepl("qPCR|GPL18522", ok$accession)], take, cap = NULL))
  rna <- comb(lapply(ok$accession[ok$parent == "GSE47774"], take, cap = cap))
  list(arr = mk(arr$X, arr$y, "array"), rna = mk(rna$X, rna$y, "rnaseq"))
}

build_metsim <- function() {
  genes <- load_index("validation")
  prep <- function(acc) {
    e <- load_expr(acc)
    ph <- readRDS(file.path(PROC_DIR, paste0(acc, "_pheno.rds")))
    i <- match(colnames(e), rownames(ph)); stopifnot(!anyNA(i)); ph <- ph[i, , drop = FALSE]
    stopifnot(all(genes %in% rownames(e)))
    bmi <- as.numeric(ph$`bmi:ch1`)
    y <- ifelse(bmi >= 30, 1L, ifelse(bmi < 25, 0L, NA_integer_))
    keep <- which(!is.na(y))
    list(X = t(rank_transform(e[genes, keep, drop = FALSE])), y = y[keep],
         id = as.character(ph$title[keep]))
  }
  list(rna = prep("GSE135134"), arr = prep("GSE70353"))
}
drop_overlap <- function(te, tr) {
  k <- !(te$id %in% tr$id)
  message("  removed ", sum(!k), " test individuals that are also in the training cohort")
  list(X = te$X[k, , drop = FALSE], y = te$y[k], id = te$id[k])
}

build_fever <- function() {
  ph <- readRDS(file.path(PROC_DIR, "GSE211567_pheno.rds"))
  e <- load_expr("GSE211567"); genes <- load_index("stress")
  stopifnot(all(colnames(e) == rownames(ph)), all(genes %in% rownames(e)))
  lab <- ph$`infection type:ch1`; site <- ph$`site:ch1`
  keep <- which(ph$geo_matched & lab %in% c("Bacterial", "Viral"))
  X <- t(rank_transform(e[genes, keep, drop = FALSE]))
  list(X = X, y = as.integer(lab[keep] == "Bacterial"), site = site[keep])
}

# ---- unit tests ------------------------------------------------------------
run_unit_tests <- function() {
  m <- matrix(c(5, 1, 3, 3, 2, 2, 2, 2), 4, 2)
  r <- rank_transform(m)
  stopifnot(all(abs(r[, 1] - c(4, 1, 2.5, 2.5) / 4) < 1e-12),
            all(abs(r[, 2] - 2.5 / 4) < 1e-12), all(r > 0 & r <= 1))
  set.seed(1); n <- 60; g <- 30
  X <- matrix(runif(n * g), n, g, dimnames = list(NULL, paste0("g", 1:g)))
  y <- rep(0:1, each = n / 2)
  X[y == 1, "g1"] <- runif(n / 2, 0, 0.4); X[y == 1, "g2"] <- runif(n / 2, 0.6, 1)
  X[y == 0, "g1"] <- runif(n / 2, 0.6, 1); X[y == 0, "g2"] <- runif(n / 2, 0, 0.4)
  fit <- ktsp_fit(X, y, k = 1, ncand = g)
  stopifnot(setequal(c(fit$pairs[[1]]$i, fit$pairs[[1]]$j), c("g1", "g2")))
  sc <- ktsp_predict(fit, X)
  stopifnot(all(sc[y == 1] > 0), all(sc[y == 0] < 0))
  message("Unit tests passed")
}
run_unit_tests()

# ---- assemble tasks --------------------------------------------------------
tasks <- list()
add_task <- function(name, direction, tr, te) {
  for (s in list(tr, te)) if (length(unique(s$y)) < 2) stop(name, " ", direction, ": a split has one class")
  tasks[[length(tasks) + 1]] <<- list(name = name, direction = direction, tr = tr, te = te)
}

sq <- build_seqc()
add_task("SEQC B-vs-A", "array -> RNA-seq", sq$arr, sq$rna)
add_task("SEQC B-vs-A", "RNA-seq -> array", sq$rna, sq$arr)

message("METSIM:")
ms <- build_metsim()
te <- drop_overlap(ms$arr, ms$rna)
add_task("METSIM obese-vs-lean", "RNA-seq -> array", mk(ms$rna$X, ms$rna$y, "rnaseq"), mk(te$X, te$y, "array"))
te <- drop_overlap(ms$rna, ms$arr)
add_task("METSIM obese-vs-lean", "array -> RNA-seq", mk(ms$arr$X, ms$arr$y, "array"), mk(te$X, te$y, "rnaseq"))

fv <- build_fever()
for (s0 in c("Sri_Lanka", "United_States")) {
  a <- fv$site == s0
  add_task("Fever bacterial-vs-viral", paste(s0, "->", setdiff(unique(fv$site), s0)),
           mk(fv$X[a, , drop = FALSE], fv$y[a], s0),
           mk(fv$X[!a, , drop = FALSE], fv$y[!a], setdiff(unique(fv$site), s0)))
}

# ---- run -------------------------------------------------------------------
results <- do.call(rbind, lapply(tasks, function(t) {
  message("\n== ", t$name, " | ", t$direction, " | train n=", length(t$tr$y),
          " (pos ", sum(t$tr$y), ") test n=", length(t$te$y), " (pos ", sum(t$te$y), ")")
  cbind(task = t$name, direction = t$direction,
        n_train = length(t$tr$y), n_test = length(t$te$y),
        run_methods(t$tr, t$te))
}))
rownames(results) <- NULL
results[, c("auc", "auc_lo", "auc_hi", "accuracy")] <-
  round(results[, c("auc", "auc_lo", "auc_hi", "accuracy")], 3)
write.csv(results, file.path(RES_DIR, "phase2_baselines.csv"), row.names = FALSE)
print(results[, c("task", "direction", "n_train", "n_test", "method",
                  "auc", "auc_lo", "auc_hi", "accuracy")], row.names = FALSE)
if (any(nzchar(results$note) & is.na(results$auc))) {
  message("\nFailed methods:"); print(results[is.na(results$auc), c("task", "direction", "method", "note")])
}
