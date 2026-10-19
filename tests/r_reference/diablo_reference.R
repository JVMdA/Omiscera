suppressMessages(library(mixOmics))

# Run from the repo root: Rscript tests/r_reference/diablo_reference.R
DIR <- "tests/fixtures/diablo_reference"

read_block <- function(name) {
  df <- read.csv(file.path(DIR, paste0("block_", name, ".csv")), row.names = 1, check.names = FALSE)
  as.matrix(df)
}
T_ <- read_block("transcriptomics")
P_ <- read_block("proteomics")
M_ <- read_block("metabolomics")
grp <- read.csv(file.path(DIR, "group.csv"), row.names = 1)
train_s <- read.csv(file.path(DIR, "train_samples.csv"))$sample
test_s <- read.csv(file.path(DIR, "test_samples.csv"))$sample

X <- list(transcriptomics = T_[train_s, ], proteomics = P_[train_s, ], metabolomics = M_[train_s, ])
Y <- factor(grp[train_s, 1])
cat("train Y levels:", levels(Y), "counts:", table(Y), "\n")

design <- matrix(0.1, nrow = 3, ncol = 3, dimnames = list(names(X), names(X)))
diag(design) <- 0

keepX <- list(transcriptomics = c(10, 8), proteomics = c(10, 8), metabolomics = c(10, 8))
ncomp <- 2

fit <- block.splsda(X = X, Y = Y, ncomp = ncomp, keepX = keepX, design = design, scheme = "horst", scale = TRUE,
                     init = "svd", tol = 1e-06, max.iter = 100, near.zero.var = FALSE)

cat("design actually used:\n"); print(fit$design)
cat("iterations:", fit$iter, "\n")

dir.create(file.path(DIR, "out"), showWarnings = FALSE)
out <- file.path(DIR, "out")

for (b in names(X)) {
  write.csv(data.frame(sample = rownames(fit$variates[[b]]), fit$variates[[b]]), file.path(out, paste0("variates_", b, ".csv")), row.names = FALSE)
  write.csv(data.frame(feature = rownames(fit$loadings[[b]]), fit$loadings[[b]]), file.path(out, paste0("loadings_", b, ".csv")), row.names = FALSE)
}
write.csv(data.frame(sample = rownames(fit$variates$Y), fit$variates$Y), file.path(out, "variates_Y.csv"), row.names = FALSE)

# AVE / explained variance per block, per component
ave_df <- do.call(rbind, lapply(names(fit$AVE$AVE_X), function(b) data.frame(block = b, comp = seq_along(fit$AVE$AVE_X[[b]]), AVE = fit$AVE$AVE_X[[b]])))
write.csv(ave_df, file.path(out, "ave.csv"), row.names = FALSE)

# weights (per-component block-Y correlation)
write.csv(data.frame(block = rownames(fit$weights), fit$weights), file.path(out, "weights.csv"), row.names = FALSE)

# selected features per block/component
sel <- do.call(rbind, lapply(names(X), function(b) {
  do.call(rbind, lapply(1:ncomp, function(h) {
    v <- selectVar(fit, block = b, comp = h)
    data.frame(block = b, comp = h, feature = v[[b]]$name, value = v[[b]]$value[, 1])
  }))
}))
write.csv(sel, file.path(out, "selected.csv"), row.names = FALSE)

# prediction on the held-out test set
Xtest <- list(transcriptomics = T_[test_s, , drop = FALSE], proteomics = P_[test_s, , drop = FALSE], metabolomics = M_[test_s, , drop = FALSE])
pr <- predict(fit, newdata = Xtest, dist = "all")

for (b in names(X)) {
  write.csv(data.frame(sample = rownames(pr$variates[[b]]), pr$variates[[b]]), file.path(out, paste0("predict_variates_", b, ".csv")), row.names = FALSE)
  for (d in c("max.dist", "centroids.dist", "mahalanobis.dist")) {
    write.csv(data.frame(sample = rownames(pr$class[[d]][[b]]), pr$class[[d]][[b]]), file.path(out, paste0("predict_class_", d, "_", b, ".csv")), row.names = FALSE)
  }
}
for (d in c("max.dist", "centroids.dist", "mahalanobis.dist")) {
  write.csv(data.frame(sample = rownames(pr$MajorityVote[[d]]), pr$MajorityVote[[d]]), file.path(out, paste0("predict_majorityvote_", d, ".csv")), row.names = FALSE)
  write.csv(data.frame(sample = rownames(pr$WeightedVote[[d]]), pr$WeightedVote[[d]]), file.path(out, paste0("predict_weightedvote_", d, ".csv")), row.names = FALSE)
}
write.csv(data.frame(sample = rownames(pr$AveragedPredict.class$max.dist), pr$AveragedPredict.class$max.dist), file.path(out, "predict_averagedpredict_max.csv"), row.names = FALSE)
write.csv(data.frame(sample = rownames(pr$WeightedPredict.class$max.dist), pr$WeightedPredict.class$max.dist), file.path(out, "predict_weightedpredict_max.csv"), row.names = FALSE)

cat("true test labels:", paste(grp[test_s, 1], collapse = ","), "\n")
write.csv(data.frame(sample = test_s, true = grp[test_s, 1]), file.path(out, "test_truth.csv"), row.names = FALSE)

# a one-component, single-block spls (2 blocks: just T and Y) as the simplest possible sanity check too
cat("DONE\n")
