# ---- SETUP ----
pkgs <- c("matrixStats", "ggplot2", "dplyr", "SeuratObject", "Seurat",
          "Matrix", "ggrepel", "patchwork", "scales", "WGCNA", "AnnotationDbi",
          "org.Mm.eg.db", "GO.db", "hdWGCNA", "limma", "rhdf5")
for (pkg in pkgs) library(pkg, character.only = TRUE)

theme_qc <- function(base_size = 10) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(linewidth = 0.25, colour = "grey92"),
      panel.border     = element_rect(colour = "grey65", fill = NA, linewidth = 0.4),
      strip.background = element_rect(fill = "grey96", colour = "grey65",
                                      linewidth = 0.4),
      strip.text       = element_text(face = "bold", size = rel(0.85)),
      plot.title       = element_text(face = "bold", size = rel(1.15)),
      plot.subtitle    = element_text(colour = "grey30", size = rel(0.90)),
      plot.caption     = element_text(colour = "grey45", size = rel(0.70),
                                      hjust = 0, margin = margin(t = 8)),
      plot.title.position   = "plot",
      plot.caption.position = "plot",
      legend.key       = element_blank(),
      axis.ticks       = element_line(linewidth = 0.3, colour = "grey65"))
}
theme_set(theme_qc())
theme_caption <- theme(
  plot.caption          = element_text(colour = "grey45", size = 7, hjust = 0,
                                       margin = margin(t = 8)),
  plot.caption.position = "plot")

cond_colours <- c(
  AL_young = scales::hue_pal()(1),                    # red
  AL_old   = scales::hue_pal(l = 20, c = 80)(1),      # dark red
  CR_young = scales::hue_pal()(7)[3],                 # green
  CR_old   = scales::hue_pal(l = 20, c = 80)(7)[3])   # dark green
condition_levels <- names(cond_colours)
condition_labels <- c(AL_young = "AL, young", AL_old = "AL, old",
                      CR_young = "CR, young", CR_old = "CR, old")

scale_cond_colour <- function(name = "Condition", ...)
  scale_colour_manual(values = cond_colours, limits = condition_levels,
                      labels = condition_labels, name = name, ...)
scale_cond_fill <- function(name = "Condition", ...)
  scale_fill_manual(values = cond_colours, limits = condition_levels,
                    labels = condition_labels, name = name, ...)

# ---- PARAMETERS ----
data_path <- "../../../../DATASETS/Mouse/Kidney/Transcriptomics/Single_cell/In_house/Diet_driven_differential_aging_AG_Ulrich_Mueller/data_depthnorm.rds"
data_path_prenorm <- "../../../../DATASETS/Mouse/Kidney/Transcriptomics/Single_cell/In_house/Diet_driven_differential_aging_AG_Ulrich_Mueller/data.rds"
mol_root <- "/cellfile/datapublic/flopes/PhD_Mouse_Kidney_snRNA_Diets_Age_only/Datasets"
out_dir  <- "../Output/00_pipeline_after_depthnorm"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

base_celltypes <- c("PT-S1", "PT-S2", "PT-S3")

n_hvg_per_sample       <- 3000   # top-N HVGs per sample before the voting
min_detection_fraction <- 0.1    # gene seen in >=10% of the cells of a sample

k_metacell       <- 5                     # cells per metacell
metacell_methods <- c("knn", "random")    # "random" as negative control
method_labels    <- c(knn    = "kNN metacells",
                      random = "random assignment metacells")
method_labeller  <- ggplot2::as_labeller(method_labels)
reference_method <- "knn"

# -metacells = built per segment then pooled, -singlecells = cells pooled then metacells built
pool_prefix_mc   <- "PT-pool-metacells"
pool_prefix_sc   <- "PT-pool-singlecells"
pool_variants    <- c("raw", "z")
pooled_celltypes <- c(paste0(pool_prefix_mc, "-", pool_variants),
                      paste0(pool_prefix_sc, "-", pool_variants))
all_celltypes    <- c(base_celltypes, pooled_celltypes)   # canonical plot order

go_min_term_size <- 15    # GO:BP term size on the org.Mm.eg.db scale
go_max_term_size <- 300   # keeps 4 of 5 PT anchor terms, drops generic parents
n_perm           <- 10000 # gene-label permutations per test

seed <- 12345
set.seed(seed)

# ---- FUNCTIONS ----
# generic helpers
with_seed <- function(seed, expr) {
  old <- if (exists(".Random.seed", .GlobalEnv)) get(".Random.seed", .GlobalEnv)
  set.seed(seed)
  on.exit(if (!is.null(old)) assign(".Random.seed", old, .GlobalEnv), add = TRUE)
  force(expr)
}
fig_registry <- list()   # every saved figure, re-rendered in THESIS FIGURES
save_pdf <- function(plots, file, width, height, encoding = "default") {
  plots <- Filter(Negate(is.null),
                  if (is.list(plots) && !inherits(plots, "gg")) plots else list(plots))
  grDevices::pdf(file.path(out_dir, file), width = width, height = height,
                 onefile = TRUE, encoding = encoding)
  on.exit(grDevices::dev.off(), add = TRUE)
  for (p in plots) print(p)
  invisible(NULL)
}
fig_caption <- function(...) {
  paste0(paste(c(...), collapse = "  |  "),
         sprintf("\nk = %d cells per metacell  |  seed = %d  |  generated %s",
                 k_metacell, seed, format(Sys.Date())))
}
load_expr_mat <- function(obj, gene_set = NULL) {
  e <- as.matrix(SeuratObject::LayerData(obj, assay = "RNA", layer = "data"))
  if (is.null(gene_set)) return(e)
  e[intersect(gene_set, rownames(e)), , drop = FALSE]
}
row_z <- function(e) {
  s <- matrixStats::rowSds(e)
  s[!is.finite(s) | s == 0] <- 1
  (e - rowMeans(e)) / s
}
apply_per_dataset <- function(nested, fn) {
  out <- do.call(rbind, lapply(seq_len(nrow(datasets)), function(i) {
    d   <- datasets[i, ]
    res <- fn(nested[[d$method]][[d$celltype]][[d$sample]], d$method, d$celltype,
              d$sample)
    if (is.null(res)) return(NULL)
    data.frame(d, res, row.names = NULL)
  }))
  out$celltype <- factor(out$celltype, levels = all_celltypes)
  out$method   <- factor(out$method,   levels = metacell_methods)
  out
}
as_metacell_object <- function(data_mat, smp) {
  placeholder <- data_mat
  placeholder[placeholder < 0] <- 0
  obj <- CreateSeuratObject(counts = as(as.matrix(placeholder), "CsparseMatrix"),
                            min.cells = 0, min.features = 0)
  SeuratObject::LayerData(obj, assay = "RNA", layer = "data") <-
    as(as.matrix(data_mat), "CsparseMatrix")
  obj$sample    <- smp
  obj$condition <- sample_to_condition[[smp]]
  obj$batch     <- sample_to_batch[[smp]]
  obj
}
metacell_groups <- function(obj, features, i, method) {
  n <- ncol(obj)
  if (method == "random")
    return(with_seed(seed + i, split(
      sample.int(n, (n %/% k_metacell) * k_metacell),
      rep(seq_len(n %/% k_metacell), each = k_metacell))))
  obj   <- ScaleData(obj, features = features, verbose = FALSE)
  n_pcs <- min(5L, n - 1L, length(features) - 1L)
  obj   <- RunPCA(obj, features = features, npcs = n_pcs, verbose = FALSE)
  obj$.hdwgcna_group <- "all"
  Idents(obj) <- ".hdwgcna_group"
  obj <- SetupForWGCNA(obj, gene_select = "custom", features = features,
                       wgcna_name = "mc")
  obj <- with_seed(seed + i, MetacellsByGroups(
    seurat_obj = obj, group.by = ".hdwgcna_group", ident.group = ".hdwgcna_group",
    k = k_metacell, min_cells = 50L, max_shared = 0L, reduction = "pca",
    dims = seq_len(n_pcs), slot = "counts", assay = "RNA", wgcna_name = "mc"))
  lapply(strsplit(as.character(GetMetacellObject(obj, wgcna_name = "mc")$cells_merged),
                  ",", fixed = TRUE),
         function(cl) match(trimws(cl), colnames(obj)))
}

# segment composition of every metacell
add_composition <- function(mc, groups, segs) {
  comp <- t(vapply(groups, function(j)
    as.numeric(table(factor(segs[j], levels = base_celltypes))) / length(j),
    numeric(length(base_celltypes))))
  dimnames(comp) <- list(colnames(mc), paste0("frac_", base_celltypes))
  mc@meta.data <- cbind(mc@meta.data, comp)
  mc
}

# single cells -> metacells: sum raw counts, log-normalise on the full gene space
build_metacells <- function(obj, gene_set, i, method, smp, segs = NULL) {
  groups <- metacell_groups(obj, intersect(gene_set, rownames(obj)), i, method)
  cnt <- SeuratObject::LayerData(obj, assay = "RNA", layer = "counts")
  mat <- vapply(groups, function(j) Matrix::rowSums(cnt[, j, drop = FALSE]),
                numeric(nrow(cnt)))
  dimnames(mat) <- list(rownames(cnt), paste0("MC_", seq_along(groups)))
  mc <- CreateSeuratObject(counts = as(mat, "CsparseMatrix"),
                           min.cells = 0, min.features = 0)
  mc <- NormalizeData(mc, verbose = FALSE)
  mc$lib_size <- Matrix::colSums(mat)
  mc <- mc[intersect(gene_set, rownames(mc)), ]
  mc$sample <- smp; mc$condition <- sample_to_condition[[smp]]
  mc$batch  <- sample_to_batch[[smp]]
  mc@misc$n_possible <- ncol(obj) %/% k_metacell
  if (is.null(segs)) mc else add_composition(mc, groups, segs)
}

# metacells on an already scaled matrix: grouping and aggregation on the same values
build_metacells_z <- function(X, i, method, smp) {
  groups <- metacell_groups(as_metacell_object(X, smp), rownames(X), i, method)
  mat <- vapply(groups, function(j) rowMeans(X[, j, drop = FALSE]), numeric(nrow(X)))
  dimnames(mat) <- list(rownames(X), paste0("MC_", seq_along(groups)))
  mc <- as_metacell_object(mat, smp)
  mc@misc$n_possible <- ncol(X) %/% k_metacell
  add_composition(mc, groups, sub("__.*", "", colnames(X)))
}

# all PT cells of one sample -> metacells
build_pool_sc <- function(gene_set, smp, i, method) {
  o <- subset(seurat_obj, cells = colnames(seurat_obj)[seurat_obj$sample == smp])
  build_metacells(o, gene_set, i, method, smp,
                  segs = as.character(o$Annotation_lvl1))
}

# gene x gene Pearson correlation across the metacells of one sample
coexpression_matrix <- function(mc, gene_set) {
  e   <- load_expr_mat(mc, gene_set)
  sds <- matrixStats::rowSds(e)
  cm  <- stats::cor(t(e[is.finite(sds) & sds > 1e-8, , drop = FALSE]),
                    method = "pearson")
  cm[is.na(cm)] <- 0
  cm
}

# cor between the two upper-triangle correlation vectors, on the common gene set
coexpression_correlation <- function(cor_a, cor_b, genes, min_genes = 50L) {
  g <- Reduce(intersect, list(genes, rownames(cor_a), rownames(cor_b)))
  if (length(g) < min_genes) return(list(r = NA_real_, n = length(g)))
  ut <- upper.tri(matrix(0, length(g), length(g)))
  list(r = stats::cor(cor_a[g, g][ut], cor_b[g, g][ut], method = "pearson"),
       n = length(g))
}

# mean |atanh(r)| of a gene to all others, diagonal excluded
connectivity_matrix <- function(cms) {
  g <- Reduce(intersect, lapply(cms, rownames))
  m <- vapply(cms, function(cm) {
    a <- abs(cm[g, g, drop = FALSE])
    diag(a) <- 0
    rowSums(atanh(pmin(a, 0.999))) / (length(g) - 1)   # 0.999 guards atanh(1) = Inf
  }, numeric(length(g)))
  dimnames(m) <- list(g, names(cms))
  m
}

# Shannon entropy H of a gene's |r| profile, top_k to a gene's k strongest partners per sample
entropy_matrix <- function(cms, top_k = NULL, keep = NULL) {
  g <- Reduce(intersect, lapply(cms, rownames))
  m <- vapply(cms, function(cm) {
    a <- abs(cm[g, g, drop = FALSE])
    diag(a) <- 0
    if (!is.null(keep)) a[!keep] <- 0
    else if (!is.null(top_k))
      a[a < matrixStats::rowOrderStats(a, which = length(g) - top_k + 1L)] <- 0
    p  <- a / rowSums(a)
    lp <- log(p)
    lp[!is.finite(lp)] <- 0              # p = 0 contributes 0, NaN stays NaN
    -rowSums(p * lp)
  }, numeric(length(g)))
  dimnames(m) <- list(g, names(cms))
  m
}

# sample-level PCA on mean expression or on the correlation vector
sample_pca <- function(objs, type, gene_set = NULL) {
  if (type == "expression") {
    g   <- Reduce(intersect, c(list(gene_set), lapply(objs, rownames)))
    mat <- t(vapply(objs, function(o) rowMeans(load_expr_mat(o, g)),
                    numeric(length(g))))
  } else {
    common <- Reduce(intersect, lapply(objs, rownames))
    ut  <- upper.tri(matrix(0, length(common), length(common)))
    mat <- t(vapply(objs, function(cm) cm[common, common][ut], numeric(sum(ut))))
  }
  sds <- matrixStats::colSds(mat)
  stats::prcomp(mat[, is.finite(sds) & sds > 1e-8, drop = FALSE],
                scale. = type == "expression")
}

# mean between-group over mean within-group distance
separation_ratio <- function(coords, grouping) {
  d    <- as.matrix(stats::dist(coords))
  same <- outer(grouping, grouping, "==")
  diag(same) <- NA
  mean(d[!same], na.rm = TRUE) / mean(d[same], na.rm = TRUE)
}

# shared-term matrix over one gene set, from the already filtered go_term2gene
go_shared_matrix <- function(gene_symbols, max_size = go_max_term_size) {
  gene_symbols <- sort(unique(gene_symbols))
  ann <- go_term2gene[
    go_term2gene$GO %in% names(go_term_size)[go_term_size <= max_size] &
      go_term2gene$SYMBOL %in% gene_symbols, ]
  term_ids  <- unique(ann$GO)
  incidence <- Matrix::sparseMatrix(
    i = match(ann$SYMBOL, gene_symbols), j = match(ann$GO, term_ids), x = 1,
    dims = c(length(gene_symbols), length(term_ids)),
    dimnames = list(gene_symbols, term_ids))
  shared_any <- as.matrix(incidence %*% Matrix::t(incidence)) > 0
  diag(shared_any) <- FALSE
  list(shared_any = shared_any, annotated = Matrix::rowSums(incidence) > 0)
}

# |cor| per gene pair, shared-term label, symmetric shared matrix
go_pair_vectors <- function(cor_mat, go) {
  genes <- intersect(rownames(cor_mat), rownames(go$shared_any))
  genes <- genes[go$annotated[genes]]
  cm <- cor_mat[genes, genes, drop = FALSE]
  sh <- go$shared_any[genes, genes, drop = FALSE]
  ut <- upper.tri(cm)
  list(genes = genes, shared_mat = sh, upper = ut,
       abs_cor = abs(as.numeric(cm[ut])), shared = as.logical(sh[ut]))
}

# enrichment of shared-GO pairs in the top-q |cor| tail, with permutation p
go_edge_enrichment <- function(cor_mat, go, top_q, n_perm, seed) {
  pv       <- go_pair_vectors(cor_mat, go)
  n_pairs  <- length(pv$abs_cor)
  n_shared <- sum(pv$shared)
  high     <- pv$abs_cor >= stats::quantile(pv$abs_cor, 1 - top_q, na.rm = TRUE)
  n_high   <- sum(high)
  
  # 2x2 table: high/low |cor| against shared/not-shared GO term
  a <- sum(high & pv$shared); b <- n_high - a
  cc <- n_shared - a;         d <- n_pairs - n_high - cc
  
  # permute gene labels, count how many of the FIXED high edges land on a shared
  # pair -- cost is O(n_high), not O(n_pairs)
  high_idx <- which(pv$upper, arr.ind = TRUE)[high, , drop = FALSE]
  ng       <- length(pv$genes)
  a_perm <- with_seed(seed, vapply(seq_len(n_perm), function(i) {
    p <- sample.int(ng)
    sum(pv$shared_mat[cbind(p[high_idx[, 1]], p[high_idx[, 2]])])
  }, numeric(1)))
  
  data.frame(a = a, b = b, c = cc, d = d, n_high = n_high, n_shared = n_shared,
             odds_ratio = (as.double(a) * d) / (as.double(b) * cc),
             perm_p = (sum(a_perm >= a) + 1) / (n_perm + 1))
}

# Mann-Whitney shift over the whole distribution instead of the tail only
go_distribution_shift <- function(cor_mat, go, n_perm, seed) {
  pv      <- go_pair_vectors(cor_mat, go)
  n_pairs <- length(pv$abs_cor)
  n_sh    <- sum(pv$shared)
  ranks   <- rank(pv$abs_cor)
  auc_of  <- function(rank_sum, s) {
    s <- as.double(s)
    (rank_sum - s * (s + 1) / 2) / (s * (n_pairs - s))
  }
  auc_obs <- auc_of(sum(ranks[pv$shared]), n_sh)
  
  # the permutation leaves n_shared fixed, so only the rank sum changes
  ng <- length(pv$genes)
  rank_mat <- matrix(0, ng, ng)
  rank_mat[pv$upper] <- ranks
  rank_mat <- rank_mat + t(rank_mat)
  shared_idx <- which(pv$upper & pv$shared_mat, arr.ind = TRUE)
  auc_perm <- with_seed(seed, vapply(seq_len(n_perm), function(i) {
    q <- sample.int(ng)     # the inverse of a uniform permutation is uniform
    auc_of(sum(rank_mat[cbind(q[shared_idx[, 1]], q[shared_idx[, 2]])]), n_sh)
  }, numeric(1)))
  
  data.frame(auc = auc_obs, perm_p = (sum(auc_perm >= auc_obs) + 1) / (n_perm + 1))
}

# odds ratio across the whole edge-fraction grid, no permutation p
go_edge_enrichment_sweep <- function(cor_mat, go, q_grid) {
  pv         <- go_pair_vectors(cor_mat, go)
  n_pairs    <- length(pv$abs_cor)
  n_shared   <- sum(pv$shared)
  ord        <- order(pv$abs_cor, decreasing = TRUE)
  cum_shared <- cumsum(pv$shared[ord])          # a as a function of n_high
  n_high <- n_pairs - findInterval(stats::quantile(pv$abs_cor, 1 - q_grid,
                                                   na.rm = TRUE),
                                   rev(pv$abs_cor[ord]), left.open = TRUE)
  a  <- cum_shared[pmax(n_high, 1L)]
  b  <- n_high - a
  cc <- n_shared - a
  d  <- n_pairs - n_high - cc
  data.frame(top_q = q_grid, odds_ratio = (as.double(a) * d) / (as.double(b) * cc))
}

# samples significantly above the null value
significance_counts <- function(df, stat, null_value, alpha = 0.05) {
  out <- dplyr::summarise(dplyr::group_by(df, method, celltype),
                          n_total = dplyr::n(),
                          n_sig = sum(perm_p < alpha & .data[[stat]] > null_value),
                          .groups = "drop")
  out$label <- paste0(out$n_sig, "/", out$n_total)
  out
}

# ---- DATA PREPARATION ----
seurat_obj <- SeuratObject::UpdateSeuratObject(readRDS(data_path))
# the metadata label the restricted diet "DR"; the thesis uses CR (caloric restriction)
seurat_obj$condition <- sub("^DR_", "CR_", as.character(seurat_obj$condition))

# library size of every cell of every annotated cell type, before any subsetting
lib_all_df <- data.frame(
  cell     = colnames(seurat_obj),
  sample   = as.character(seurat_obj$sample),
  celltype = as.character(seurat_obj$Annotation_lvl1),
  lib      = Matrix::colSums(SeuratObject::LayerData(seurat_obj, assay = "RNA",
                                                     layer = "counts")),
  row.names = NULL, stringsAsFactors = FALSE)

# gene universe over every annotated cell type, for the section A GO panel
genes_kidney_expressed <- rownames(seurat_obj)[Matrix::rowSums(
  SeuratObject::LayerData(seurat_obj, assay = "RNA", layer = "counts")) > 0]

seurat_obj <- subset(seurat_obj, cells = colnames(seurat_obj)[
  seurat_obj$Annotation_lvl1 %in% base_celltypes])

sample_meta   <- unique(seurat_obj@meta.data[, c("sample", "condition", "batch")])
sample_meta[] <- lapply(sample_meta, as.character)
sample_meta   <- sample_meta[order(factor(sample_meta$condition,
                                          levels = condition_levels),
                                   sample_meta$sample), ]
rownames(sample_meta) <- NULL
sample_to_condition <- setNames(sample_meta$condition, sample_meta$sample)
sample_to_batch     <- setNames(sample_meta$batch,     sample_meta$sample)
sample_levels       <- sample_meta$sample
batch_levels        <- sort(unique(sample_meta$batch))
batch_shapes        <- setNames(c(16, 17, 15, 18, 8, 4, 3, 7)[seq_along(batch_levels)],
                                batch_levels)

# cell count equalisation within each cell type
counts_before <- as.data.frame(table(celltype = as.character(seurat_obj$Annotation_lvl1),
                                     sample   = as.character(seurat_obj$sample)),
                               responseName = "count")

cells_kept <- with_seed(seed, unlist(lapply(base_celltypes, function(ct) {
  in_ct  <- seurat_obj$Annotation_lvl1 == ct
  by_smp <- split(colnames(seurat_obj)[in_ct], as.character(seurat_obj$sample[in_ct]))
  unlist(lapply(by_smp, function(cells) sample(cells, min(lengths(by_smp)))))
})))
seurat_obj <- subset(seurat_obj, cells = unname(cells_kept))

counts_after <- as.data.frame(table(celltype = as.character(seurat_obj$Annotation_lvl1),
                                    sample   = as.character(seurat_obj$sample)),
                              responseName = "count")

# same quantity before the depth normalisation, cells matched by barcode
lib_prenorm_df <- local({
  o  <- SeuratObject::UpdateSeuratObject(readRDS(data_path_prenorm))
  df <- data.frame(
    cell   = colnames(o),
    sample = as.character(o$sample),
    lib    = Matrix::colSums(SeuratObject::LayerData(o, assay = "RNA",
                                                     layer = "counts")),
    row.names = NULL, stringsAsFactors = FALSE)
  rm(o); gc()
  df
})

# split into cell type x sample
objs <- setNames(lapply(base_celltypes, function(ct)
  setNames(lapply(sample_levels, function(s)
    subset(seurat_obj, cells = colnames(seurat_obj)[
      seurat_obj$Annotation_lvl1 == ct & seurat_obj$sample == s])),
    sample_levels)), base_celltypes)

# molecule_info.h5 of the libraries before the depth normalisation
mol_meta <- data.frame(
  sid    = c("SID115156", "SID115157", "SID118986", "SID118987",
             "SID136207", "SID136208", "SID136210", "SID136209"),
  sample = c("01", "02", "03", "04", "05", "06", "07", "08"),
  stringsAsFactors = FALSE)
mol_meta <- mol_meta[match(sample_levels, mol_meta$sample), ]

mol_found <- Sys.glob(file.path(mol_root, "*", "*", "cellranger", "*", "outs",
                                "molecule_info.h5"))
mol_sid   <- vapply(strsplit(sub(paste0("^", mol_root, "/"), "", mol_found), "/",
                             fixed = TRUE), `[`, character(1), 4)
mol_paths <- setNames(mol_found[match(mol_meta$sid, mol_sid)], mol_meta$sample)

# reads per molecule, restricted to the cell-associated barcodes
mol_stats <- as.data.frame(t(vapply(mol_paths, function(path) {
  h <- rhdf5::h5read(path, "/")
  on.exit(rhdf5::h5closeAll(), add = TRUE)
  pf <- h$barcode_info$pass_filter              # (barcode_idx, library_idx, genome_idx)
  if (nrow(pf) == 3L && ncol(pf) != 3L) pf <- t(pf)
  cells <- unique(as.vector(pf[, 1]))
  r <- as.numeric(h$count[h$barcode_idx %in% cells])
  c(reads_pc = sum(r) / length(cells), singleton = mean(r == 1))
}, numeric(2))))
mol_stats$sample    <- factor(rownames(mol_stats), levels = sample_levels)
mol_stats$condition <- factor(sample_to_condition[rownames(mol_stats)],
                              levels = condition_levels)
rownames(mol_stats) <- NULL
mol_target <- min(mol_stats$reads_pc)   # aggr equalises to the shallowest library

# ---- GENE SELECTION ----
# detection floor, then top-N HVGs per sample, both voted within condition
select_genes_per_celltype <- function(objs, n_hvg, min_detection) {
  sample_ids <- names(objs)
  conditions <- unique(sample_to_condition[sample_ids])
  samples_of <- function(cond) sample_ids[sample_to_condition[sample_ids] == cond]
  
  detection <- do.call(cbind, lapply(objs, function(o)
    Matrix::rowMeans(SeuratObject::LayerData(o, assay = "RNA", layer = "data") > 0)))
  colnames(detection) <- sample_ids
  detected <- unique(unlist(lapply(conditions, function(cond) {
    sub <- detection[, samples_of(cond), drop = FALSE]
    rownames(sub)[rowSums(sub >= min_detection) == ncol(sub)]
  })))
  
  hvg_info <- lapply(objs, function(o) {
    cnt <- SeuratObject::LayerData(o, assay = "RNA", layer = "counts")
    tmp <- CreateSeuratObject(counts = cnt[intersect(detected, rownames(cnt)), ,
                                           drop = FALSE],
                              min.cells = 0, min.features = 0)
    tmp <- FindVariableFeatures(tmp, selection.method = "vst", nfeatures = n_hvg,
                                verbose = FALSE)
    hv <- SeuratObject::HVFInfo(tmp, assay = "RNA", method = "vst")
    names(hv) <- sub(".*_", "", names(hv))      # v5 layer prefixes, if present
    list(hvg = VariableFeatures(tmp),
         vst = hv[order(hv$variance.standardized, decreasing = TRUE), , drop = FALSE])
  })
  
  list(selected = unique(unlist(lapply(conditions, function(cond)
    Reduce(intersect, lapply(hvg_info[samples_of(cond)], `[[`, "hvg"))))),
    detected  = detected,
    expressed = rownames(detection)[rowSums(detection > 0) > 0],
    vst       = lapply(hvg_info, `[[`, "vst"))
}

gene_selection <- setNames(lapply(base_celltypes, function(ct)
  select_genes_per_celltype(objs[[ct]], n_hvg_per_sample, min_detection_fraction)),
  base_celltypes)

pt_gene_set <- Reduce(intersect, lapply(gene_selection, `[[`, "selected"))
gene_sets   <- setNames(rep(list(pt_gene_set), length(all_celltypes)), all_celltypes)

# ---- METACELL CONSTRUCTION ----
datasets <- expand.grid(sample = sample_levels, celltype = all_celltypes,
                        method = metacell_methods, stringsAsFactors = FALSE)
datasets$condition <- factor(sample_to_condition[datasets$sample],
                             levels = condition_levels)

metacells <- setNames(vector("list", length(metacell_methods)), metacell_methods)

for (method in metacell_methods) {
  metacells[[method]] <- setNames(lapply(base_celltypes, function(ct)
    setNames(lapply(seq_along(sample_levels), function(i)
      build_metacells(objs[[ct]][[i]], gene_sets[[ct]], i, method, sample_levels[i])),
      sample_levels)), base_celltypes)
  
  n_poss <- function(smp) sum(vapply(base_celltypes, function(ct)
    metacells[[method]][[ct]][[smp]]@misc$n_possible, numeric(1)))
  
  # pool-metacells-raw: segment metacell counts concatenated, normalised after
  metacells[[method]][[paste0(pool_prefix_mc, "-raw")]] <-
    setNames(lapply(sample_levels, function(smp) {
      mats <- lapply(base_celltypes, function(ct) {
        cnt <- as.matrix(SeuratObject::LayerData(metacells[[method]][[ct]][[smp]],
                                                 assay = "RNA", layer = "counts"))
        colnames(cnt) <- paste0(ct, "__", colnames(cnt))
        cnt
      })
      common <- Reduce(intersect, lapply(mats, rownames))
      X  <- do.call(cbind, lapply(mats, function(m) m[common, , drop = FALSE]))
      mc <- CreateSeuratObject(counts = as(X, "CsparseMatrix"),
                               min.cells = 0, min.features = 0)
      mc <- NormalizeData(mc, verbose = FALSE)         # denominator = full gene set
      mc$lib_size <- Matrix::colSums(X)
      mc <- mc[intersect(pt_gene_set, rownames(mc)), ] # subset only afterwards
      mc$sample <- smp; mc$condition <- sample_to_condition[[smp]]
      mc$batch  <- sample_to_batch[[smp]]
      mc@misc$n_possible <- n_poss(smp)
      mc
    }), sample_levels)
  
  # pool-metacells-z: z-scored per segment on the data layer, then concatenated
  metacells[[method]][[paste0(pool_prefix_mc, "-z")]] <-
    setNames(lapply(sample_levels, function(smp) {
      mats <- lapply(base_celltypes, function(ct) {
        e <- row_z(load_expr_mat(metacells[[method]][[ct]][[smp]], pt_gene_set))
        colnames(e) <- paste0(ct, "__", colnames(e))
        e
      })
      common <- Reduce(intersect, lapply(mats, rownames))
      mc <- as_metacell_object(
        do.call(cbind, lapply(mats, function(m) m[common, , drop = FALSE])), smp)
      mc@misc$n_possible <- n_poss(smp)
      mc
    }), sample_levels)
  
  # pool-singlecells-raw: all segments of one sample in one object, metacells on counts
  metacells[[method]][[paste0(pool_prefix_sc, "-raw")]] <-
    setNames(lapply(seq_along(sample_levels), function(i)
      build_pool_sc(pt_gene_set, sample_levels[i], i, method)), sample_levels)
  
  # pool-singlecells-z: z-scored per segment before the grouping
  metacells[[method]][[paste0(pool_prefix_sc, "-z")]] <-
    setNames(lapply(seq_along(sample_levels), function(i) {
      mats <- lapply(base_celltypes, function(ct) {
        e <- row_z(load_expr_mat(objs[[ct]][[sample_levels[i]]], pt_gene_set))
        colnames(e) <- paste0(ct, "__", colnames(e))
        e
      })
      common <- Reduce(intersect, lapply(mats, rownames))
      build_metacells_z(do.call(cbind, lapply(mats, function(m)
        m[common, , drop = FALSE])), i, method, sample_levels[i])
    }), sample_levels)
}

# equalise the metacell number across samples and methods; the four pooled arms
# additionally share one minimum so that they stay comparable to each other
n_eq_per_celltype <- vapply(all_celltypes, function(ct)
  min(vapply(metacell_methods, function(m)
    min(vapply(metacells[[m]][[ct]], ncol, numeric(1))), numeric(1))), numeric(1))
n_eq_per_celltype[pooled_celltypes] <- min(n_eq_per_celltype[pooled_celltypes])

for (method in metacell_methods) for (ct in all_celltypes) {
  metacells[[method]][[ct]] <- setNames(lapply(
    seq_along(metacells[[method]][[ct]]), function(i) {
      mc   <- metacells[[method]][[ct]][[i]]
      keep <- with_seed(seed + i, sort(sample.int(ncol(mc), n_eq_per_celltype[[ct]])))
      out  <- mc[, colnames(mc)[keep]]
      out@misc <- mc@misc                       # subset drops misc
      out@misc$n_built <- ncol(mc)
      out
    }), sample_levels)
}

# ---- CO-EXPRESSION MATRICES ----
cor_mats <- setNames(lapply(metacell_methods, function(method)
  setNames(lapply(all_celltypes, function(ct)
    lapply(metacells[[method]][[ct]], coexpression_matrix,
           gene_set = gene_sets[[ct]])), all_celltypes)), metacell_methods)

# ---- GO ANNOTATION ----
# GO:BP term -> gene table on the full mouse annotation, ancestor-propagated
go_ontology <- suppressMessages(AnnotationDbi::select(
  GO.db, keys = AnnotationDbi::keys(GO.db, keytype = "GOID"),
  columns = "ONTOLOGY", keytype = "GOID"))
go_terms <- intersect(
  go_ontology$GOID[!is.na(go_ontology$ONTOLOGY) & go_ontology$ONTOLOGY == "BP"],
  AnnotationDbi::mappedkeys(org.Mm.egGO2ALLEGS))

go_entrez <- lapply(AnnotationDbi::mget(go_terms, org.Mm.egGO2ALLEGS), unique)
go_entrez <- go_entrez[lengths(go_entrez) >= go_min_term_size &
                         lengths(go_entrez) <= go_max_term_size]
go_term_size <- lengths(go_entrez)

go_flat   <- unlist(go_entrez, use.names = FALSE)
go_symbol <- suppressMessages(AnnotationDbi::mapIds(
  org.Mm.eg.db, unique(go_flat), "SYMBOL", "ENTREZID"))

go_term2gene <- unique(data.frame(
  GO = rep(names(go_entrez), lengths(go_entrez)),
  SYMBOL = unname(go_symbol[go_flat]), stringsAsFactors = FALSE))
go_term2gene <- go_term2gene[!is.na(go_term2gene$SYMBOL), ]

go_term2name <- suppressMessages(AnnotationDbi::select(
  GO.db, keys = unique(go_term2gene$GO), columns = "TERM", keytype = "GOID"))
names(go_term2name) <- c("GO", "TERM")

# ---- A INPUT DESCRIPTION ----
cell_count_df <- merge(counts_before, counts_after, by = c("celltype", "sample"),
                       suffixes = c("_before", "_after"))
cell_count_df$sample    <- factor(as.character(cell_count_df$sample),
                                  levels = sample_levels)
cell_count_df$celltype  <- factor(as.character(cell_count_df$celltype),
                                  levels = base_celltypes)
cell_count_df$condition <- factor(sample_to_condition[as.character(cell_count_df$sample)],
                                  levels = condition_levels)

p_cells <- ggplot(cell_count_df, aes(x = sample)) +
  geom_col(aes(y = count_before), fill = "grey88", width = 0.75) +
  geom_col(aes(y = count_after, fill = condition), width = 0.75) +
  geom_text(aes(y = count_after, label = count_after), vjust = 1.3, size = 2.3,
            colour = "white") +
  geom_text(aes(y = count_before,
                label = sprintf("%.0f%%", 100 * count_after / count_before)),
            vjust = -0.5, size = 2, colour = "grey45") +
  facet_wrap(~ celltype, nrow = 1) +
  scale_cond_fill() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12)), labels = scales::comma) +
  labs(x = NULL, y = "Cells",
       title = "Cell numbers equalisation across samples within each cell type",
       subtitle = "Grey = cells available before downsampling; coloured = cells kept after downsampling",
       caption = fig_caption()) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
        legend.position = "bottom")

save_pdf(p_cells, "A1_cell_counts.pdf",
         width = 3.4 * length(base_celltypes) + 1, height = 5)

# gene selection cascade
stage_labels <- c(
  expressed = "Any expressed\n(detected in >0 cells)",
  detected  = sprintf("Detected\n(>=%.0f%% of cells in every\nsample of >=1 condition)",
                      100 * min_detection_fraction),
  selected  = sprintf("Variable (HVG)\n(top-%d per sample in every\nsample of >=1 condition)",
                      n_hvg_per_sample),
  final     = "Final gene set\n(intersected across cell types)")

gene_stage_df <- do.call(rbind, lapply(base_celltypes, function(ct) {
  gs <- gene_selection[[ct]]
  n  <- c(expressed = length(gs$expressed), detected = length(gs$detected),
          selected  = length(gs$selected),  final = length(pt_gene_set))
  data.frame(celltype = ct, stage = names(n), n = as.integer(n),
             stringsAsFactors = FALSE)
}))
gene_stage_df$celltype <- factor(gene_stage_df$celltype, levels = base_celltypes)
gene_stage_df$stage    <- factor(gene_stage_df$stage, levels = names(stage_labels),
                                 labels = stage_labels)

p_genes <- ggplot(gene_stage_df, aes(x = celltype, y = n, fill = stage)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_text(aes(label = scales::comma(n)), position = position_dodge(width = 0.8),
            vjust = -0.5, size = 2.8, colour = "grey20") +
  scale_fill_manual(values = setNames(c("grey82", "#9EC9C0", "#4E8FA6", "#1F3F5B"),
                                      stage_labels), name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12)), labels = scales::comma) +
  labs(x = NULL, y = "Genes",
       title = "Gene selection pipeline",
       subtitle = "Voting within each condition allows detection of age / diet specific effects") +
  theme(legend.position   = "bottom",
        legend.text       = element_text(size = 6.5, lineheight = 1.2,
                                         margin = margin(l = 2, r = 10)),
        legend.key.height = unit(24, "pt"))

# HVG cut-off sweep
hvg_sweep_n <- c(500, 1000, 1500, 2000, 2500, 3000)   # x axis of the yield curve
hvg_run_n   <- c(1000, 1500, 2000, 2500, 3000)        # full rebuild each, expensive

hvg_curve_df <- do.call(rbind, lapply(base_celltypes, function(ct)
  do.call(rbind, lapply(names(gene_selection[[ct]]$vst), function(s) {
    v <- gene_selection[[ct]]$vst[[s]]$variance.standardized
    data.frame(celltype = ct, sample = s, rank = seq_along(v), var_std = v,
               row.names = NULL, stringsAsFactors = FALSE)
  }))))
hvg_curve_df$celltype  <- factor(hvg_curve_df$celltype, levels = base_celltypes)
hvg_curve_df$condition <- factor(sample_to_condition[hvg_curve_df$sample],
                                 levels = condition_levels)

# gene sets and yield of the full voting chain at each candidate cut-off
hvg_sets <- setNames(lapply(hvg_sweep_n, function(n) {
  sel <- lapply(setNames(base_celltypes, base_celltypes), function(ct) {
    tops  <- lapply(gene_selection[[ct]]$vst, function(v) head(rownames(v), n))
    conds <- unique(sample_to_condition[names(tops)])
    unique(unlist(lapply(conds, function(cd) Reduce(intersect,
                                                    tops[names(tops)[sample_to_condition[names(tops)] == cd]]))))
  })
  c(sel, list(final = Reduce(intersect, sel)))
}), hvg_sweep_n)

hvg_yield_df <- do.call(rbind, lapply(hvg_sweep_n, function(n)
  data.frame(n_hvg = n, celltype = c(base_celltypes, "final"),
             n_genes = unname(lengths(hvg_sets[[as.character(n)]])[
               c(base_celltypes, "final")]),
             row.names = NULL, stringsAsFactors = FALSE)))

# per-cut-off yield annotated on the topmost curve of each facet
hvg_final_n <- setNames(hvg_yield_df$n_genes[hvg_yield_df$celltype == "final"],
                        hvg_yield_df$n_hvg[hvg_yield_df$celltype == "final"])

hvg_elbow_lab <- do.call(rbind, lapply(base_celltypes, function(ct) {
  d    <- hvg_curve_df[hvg_curve_df$celltype == ct & hvg_curve_df$rank %in% hvg_sweep_n, ]
  y_hi <- tapply(d$var_std, d$rank, max)
  y_lo <- tapply(d$var_std, d$rank, min)
  n_ct <- setNames(hvg_yield_df$n_genes[hvg_yield_df$celltype == ct],
                   hvg_yield_df$n_hvg[hvg_yield_df$celltype == ct])
  data.frame(celltype = ct, rank = as.integer(names(y_hi)),
             var_std_hi = as.numeric(y_hi), var_std_lo = as.numeric(y_lo),
             label_ct    = scales::comma(n_ct[names(y_hi)]),
             label_final = scales::comma(hvg_final_n[names(y_hi)]),
             row.names = NULL, stringsAsFactors = FALSE)
}))
hvg_elbow_lab$celltype <- factor(hvg_elbow_lab$celltype, levels = base_celltypes)

p_hvg_elbow <- ggplot(hvg_curve_df, aes(rank, var_std, colour = condition,
                                        group = sample)) +
  geom_vline(xintercept = hvg_sweep_n, linetype = "dotted", colour = "grey60",
             linewidth = 0.3) +
  geom_line(linewidth = 0.4, alpha = 0.85) +
  geom_point(data = hvg_curve_df[hvg_curve_df$rank %in% hvg_sweep_n, ], size = 1.1) +
  geom_text(data = hvg_elbow_lab, aes(rank, var_std_hi, label = label_ct),
            inherit.aes = FALSE, vjust = -1.2, size = 2.1, colour = "grey20") +
  geom_text(data = hvg_elbow_lab, aes(rank, var_std_lo, label = label_final),
            inherit.aes = FALSE, vjust = 2.1, size = 2.1, colour = "grey55") +
  facet_wrap(~ celltype, nrow = 1) +
  scale_x_continuous(breaks = c(hvg_sweep_n, max(hvg_curve_df$rank))) +
  scale_y_log10(expand = expansion(mult = c(0.16, 0.14))) +
  scale_cond_colour() +
  labs(x = "Gene rank within the sample", y = "vst standardised variance (log)",
       title = "HVG cut-off: variance decay and resulting gene set size",
       subtitle = "One curve per sample, genes ordered by vst standardised variance; Dotted lines = candidate cut-offs; Upper number = selected genes within this cell type, lower number = after intersecting across cell types") +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 6, angle = 45, hjust = 1))

# full rebuild per cut-off: gene set -> kNN partition -> correlations -> connectivity
hvg_run_keys <- as.character(sort(hvg_run_n))
hvg_run_pairs <- lapply(seq_len(length(hvg_run_keys) - 1L), function(i)
  hvg_run_keys[c(i, i + 1L)])
names(hvg_run_pairs) <- vapply(hvg_run_pairs, function(p)
  sprintf("top-%s / top-%s", p[1], p[2]), character(1))

hvg_run_mc <- setNames(lapply(hvg_run_keys, function(n)
  setNames(lapply(seq_along(sample_levels), function(i)
    build_pool_sc(hvg_sets[[n]]$final, sample_levels[i], i, "knn")), sample_levels)),
  hvg_run_keys)

hvg_run_built <- vapply(hvg_run_mc, function(m) vapply(m, ncol, numeric(1)),
                        numeric(length(sample_levels)))
dimnames(hvg_run_built) <- list(sample_levels, hvg_run_keys)

# fewer metacells -> noisier correlations -> number is equalised within a pair
hvg_pair_conn <- lapply(hvg_run_pairs, function(p) {
  n_eq <- min(hvg_run_built[, p])
  out  <- lapply(setNames(p, p), function(n) connectivity_matrix(
    setNames(lapply(seq_along(sample_levels), function(i) {
      mc   <- hvg_run_mc[[n]][[i]]
      keep <- with_seed(seed + i, sort(sample.int(ncol(mc), n_eq)))
      coexpression_matrix(mc[, colnames(mc)[keep]], gene_set = hvg_sets[[n]]$final)
    }), sample_levels)))
  attr(out, "n_eq") <- n_eq
  out
})

hvg_conn_df <- do.call(rbind, lapply(names(hvg_run_pairs), function(lab) {
  p  <- hvg_run_pairs[[lab]]
  cn <- hvg_pair_conn[[lab]]
  a  <- cn[[p[1]]]; b <- cn[[p[2]]]
  g  <- intersect(rownames(a), rownames(b))
  data.frame(pair = lab, key_a = p[1], key_b = p[2],
             n_a = nrow(a), n_b = nrow(b), n_shared = length(g),
             n_eq = attr(cn, "n_eq"),
             gene = rep(g, times = length(sample_levels)),
             sample = rep(sample_levels, each = length(g)),
             conn_a = as.numeric(a[g, sample_levels]),
             conn_b = as.numeric(b[g, sample_levels]),
             row.names = NULL, stringsAsFactors = FALSE)
}))
hvg_conn_df$pair      <- factor(hvg_conn_df$pair, levels = names(hvg_run_pairs))
hvg_conn_df$condition <- factor(sample_to_condition[hvg_conn_df$sample],
                                levels = condition_levels)

hvg_conn_fit <- do.call(rbind, lapply(split(hvg_conn_df, hvg_conn_df$pair), function(d) {
  cf <- coef(lm(conn_b ~ conn_a, data = d))
  data.frame(pair = d$pair[1], slope = cf[2], intercept = cf[1],
             row.names = NULL)
}))

hvg_conn_labs <- local({
  d <- unique(hvg_conn_df[, c("pair", "key_a", "key_b", "n_a", "n_b", "n_shared", "n_eq")])
  f <- hvg_conn_fit[match(d$pair, hvg_conn_fit$pair), ]
  setNames(sprintf("x = top-%s HVGs (%s genes)  |  y = top-%s HVGs (%s genes)\n%s shared genes  |  %d metacells per sample  |  OLS slope = %.3f",
                   d$key_a, scales::comma(d$n_a), d$key_b, scales::comma(d$n_b),
                   scales::comma(d$n_shared), d$n_eq, f$slope),
           as.character(d$pair))
})
hvg_conn_lim <- range(c(hvg_conn_df$conn_a, hvg_conn_df$conn_b))

p_hvg_conn <- ggplot(hvg_conn_df, aes(conn_a, conn_b, colour = condition)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey45") +
  geom_point(size = 0.5, alpha = 0.35) +
  geom_smooth(aes(group = 1), method = "lm", formula = y ~ x, se = FALSE,
              colour = "grey55", linewidth = 0.5) +
  facet_wrap(~ pair, nrow = 1, labeller = labeller(pair = hvg_conn_labs)) +
  scale_x_continuous(limits = hvg_conn_lim) +
  scale_y_continuous(limits = hvg_conn_lim) +
  scale_cond_colour() +
  coord_fixed() +
  guides(colour = guide_legend(override.aes = list(size = 2, alpha = 1))) +
  labs(x = "Connectivity | smaller cut-off",
       y = "Connectivity | larger cut-off",
       subtitle = sprintf("One point per gene x sample on the genes shared by both sets in %s. Dashed line = identity, solid line = OLS fit pooled over samples",
                          paste0(pool_prefix_sc, "-raw"))) +
  theme(legend.position = "bottom",
        strip.text = element_text(size = 7, lineheight = 1.3))

page_w       <- 4.3 * length(base_celltypes) + 1
conn_panel_w <- (page_w - 0.95 - 0.06 * (length(hvg_run_pairs) - 1L)) /
  length(hvg_run_pairs)
row_h        <- c(5, 6, conn_panel_w + 1.5)

save_pdf(patchwork::wrap_plots(p_genes, p_hvg_elbow, p_hvg_conn,
                               ncol = 1, heights = row_h) +
           patchwork::plot_annotation(
             caption = fig_caption(
               "connectivity = mean |atanh(Pearson r)| of a gene to its neighbors"),
             theme = theme_caption),
         "A2_gene_selection.pdf", width = page_w, height = sum(row_h) + 0.2)

# GO enrichment of the selected gene set under three candidate universes
go_enrich_min_universe <- 5   # floor on the universe scale, applied by enricher()

go_universes <- list(
  kidney    = genes_kidney_expressed,
  expressed = Reduce(intersect, lapply(gene_selection, `[[`, "expressed")),
  detected  = Reduce(intersect, lapply(gene_selection, `[[`, "detected")))

frac <- function(x) vapply(strsplit(x, "/", fixed = TRUE),
                           function(p) as.numeric(p[1]) / as.numeric(p[2]), numeric(1))

enrich_df <- do.call(rbind, lapply(names(go_universes), function(u) {
  uni <- go_universes[[u]]
  d   <- as.data.frame(clusterProfiler::enricher(
    gene = intersect(pt_gene_set, uni), universe = uni,
    TERM2GENE = go_term2gene, TERM2NAME = go_term2name,
    minGSSize = go_enrich_min_universe, maxGSSize = length(uni),
    pvalueCutoff = 1, qvalueCutoff = 1, pAdjustMethod = "BH"))
  data.frame(universe = u, id = d$ID, term = d$Description,
             fold = frac(d$GeneRatio) / frac(d$BgRatio),
             count = d$Count, padj = d$p.adjust,
             n_tested = nrow(d), n_universe = length(uni),
             n_query = as.numeric(sub(".*/", "", d$GeneRatio[1])),
             stringsAsFactors = FALSE)
}))

# union of the top terms per universe, so a term missing in one panel stays visible
top_ids <- unique(unlist(lapply(split(enrich_df, enrich_df$universe),
                                function(d) head(d$id[order(d$padj)], 5))))
plot_df <- enrich_df[enrich_df$id %in% top_ids, ]

# y order follows the "detected" universe; terms absent there go to the bottom
ord   <- setNames(rep(Inf, length(top_ids)), top_ids)
d_det <- enrich_df[enrich_df$universe == "detected", ]
hit   <- intersect(top_ids, d_det$id)
ord[hit] <- d_det$padj[match(hit, d_det$id)]
id2term  <- setNames(make.unique(enrich_df$term[match(top_ids, enrich_df$id)]), top_ids)
plot_df$label <- factor(id2term[plot_df$id],
                        levels = unname(id2term[names(sort(ord, decreasing = TRUE))]))

uni_sum  <- unique(enrich_df[, c("universe", "n_universe", "n_tested", "n_query")])
n_sig    <- tapply(enrich_df$padj < 0.05, enrich_df$universe, sum)
uni_desc <- c(
  kidney    = sprintf("any expressed gene, all %d annotated cell types",
                      length(unique(lib_all_df$celltype))),
  expressed = sprintf("any expressed gene, %s only",
                      paste(base_celltypes, collapse = " / ")),
  detected  = sprintf("detected in >=%.0f%% of the %s cells of >=1 condition",
                      100 * min_detection_fraction,
                      paste(base_celltypes, collapse = " / ")))
uni_labels <- setNames(
  sprintf("Universe: %s\n%s genes  |  %s terms tested\n%d with FDR < 0.05  |  %s of %d selected genes annotated",
          uni_desc[uni_sum$universe], scales::comma(uni_sum$n_universe),
          scales::comma(uni_sum$n_tested), n_sig[uni_sum$universe],
          scales::comma(uni_sum$n_query), length(pt_gene_set)),
  uni_sum$universe)
plot_df$universe <- factor(plot_df$universe, levels = names(go_universes))

p_go_sel <- ggplot(plot_df, aes(fold, label)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey55") +
  geom_point(aes(size = count, colour = -log10(padj))) +
  facet_wrap(~ universe, ncol = length(go_universes),
             labeller = labeller(universe = uni_labels)) +
  scale_colour_gradient(low = "grey78", high = "#B2182B",
                        name = expression(-log[10]~"FDR")) +
  scale_size_continuous(range = c(1.5, 5), name = "Selected genes\nin term") +
  scale_y_discrete(labels = scales::label_wrap(40)) +
  labs(x = "Fold enrichment (gene ratio / background ratio)", y = NULL,
       title = sprintf("GO:BP enrichment of the %d selected genes under %d candidate universes",
                       length(pt_gene_set), length(go_universes)),
       subtitle = "Union of the top 5 terms of all panels, y order = p.adj in the detected universe",
       caption = fig_caption(
         sprintf("term size %d-%d on the org.Mm.eg.db scale (org.Mm.egGO2ALLEGS, ancestor-propagated), additional floor of %d genes on the universe scale",
                 go_min_term_size, go_max_term_size, go_enrich_min_universe),
         "hypergeometric test, BH adjustment")) +
  theme(legend.position = "right",
        axis.text.y = element_text(size = 7),
        strip.text  = element_text(size = 7, face = "plain", lineheight = 1.3,
                                   margin = margin(4, 3, 4, 3)))

save_pdf(p_go_sel, "A3_GO_gene_selection.pdf",
         width  = 5.4 * length(go_universes) + 4.3,
         height = 0.28 * nlevels(plot_df$label) + 3.2)

# ---- B METACELL DIAGNOSTICS ----
mc_count_df <- apply_per_dataset(metacells, function(mc, method, ct, smp)
  data.frame(n_kept = ncol(mc), n_built = mc@misc$n_built,
             n_possible = mc@misc$n_possible))

p_mc <- ggplot(mc_count_df, aes(x = sample)) +
  geom_col(aes(y = n_possible), fill = "grey88", width = 0.8) +
  geom_col(aes(y = n_built), fill = "grey68", width = 0.8) +
  geom_col(aes(y = n_kept, fill = condition), width = 0.8) +
  geom_text(aes(y = n_kept, label = n_kept), vjust = 1.3, size = 2.3,
            colour = "white") +
  geom_text(aes(y = n_built, label = sprintf("%.0f%%", 100 * n_kept / n_built)),
            vjust = -0.4, size = 1.7, colour = "grey35") +
  facet_grid(method ~ celltype, labeller = labeller(method = method_labeller)) +
  scale_cond_fill() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
  labs(x = NULL, y = "Metacells per sample",
       title = "Metacells per sample, cell type and pooling strategy",
       subtitle = "Grey = possible if every cell were used; mid grey = actually built; coloured = kept after equalisation across samples and methods") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 5.5),
        legend.position = "none")

comp_df <- apply_per_dataset(metacells, function(mc, method, ct, smp) {
  cols <- paste0("frac_", base_celltypes)
  if (!all(cols %in% colnames(mc@meta.data))) return(NULL)
  data.frame(purity = matrixStats::rowMaxs(as.matrix(mc@meta.data[, cols])))
})

# expected largest segment fraction when k cells are drawn at random
segment_frac    <- tapply(counts_after$count, counts_after$celltype, sum)[base_celltypes]
segment_frac    <- segment_frac / sum(segment_frac)
purity_comp     <- expand.grid(rep(list(0:k_metacell), length(base_celltypes)))
purity_comp     <- purity_comp[rowSums(purity_comp) == k_metacell, , drop = FALSE]
purity_prob     <- apply(purity_comp, 1, function(x)
  factorial(k_metacell) / prod(factorial(x)) * prod(segment_frac^x))
purity_baseline <- sum(purity_prob * apply(purity_comp, 1, max)) / k_metacell

p_purity <- ggplot(comp_df, aes(sample, purity, fill = condition)) +
  geom_violin(scale = "width", linewidth = 0.2, colour = "grey45") +
  geom_hline(yintercept = purity_baseline, linetype = "dashed", colour = "grey55") +
  facet_grid(method ~ celltype, labeller = labeller(method = method_labeller)) +
  scale_cond_fill() +
  scale_y_continuous(limits = c(0.25, 1)) +
  labs(x = NULL, y = "Largest segment\nfraction per metacell",
       subtitle = sprintf("Segment composition of the single-cell pools. 1 = metacell from one segment only. Dashed at %.3f = expected largest fraction when %d cells are drawn at random from a pool of %s. In the metacell pools every metacell is pure by construction.",
                          purity_baseline, k_metacell,
                          paste(sprintf("%.0f%% %s", 100 * segment_frac, base_celltypes),
                                collapse = " / "))) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 5.5),
        legend.position = "bottom")

save_pdf(patchwork::wrap_plots(p_mc, p_purity, ncol = 1, heights = c(1.4, 1)) +
           patchwork::plot_annotation(
             caption = fig_caption(
               "kNN grouping by hdWGCNA on 5 PCs of the scaled log values, max_shared = 0; ",
               "pool-metacells = metacells built per segment, then concatenated",
               "pool-singlecells = single cells concatenated, then metacells built on the pool",
               "-raw = counts summed and log-normalised after pooling",
               "-z = per segment z-scored data layer, aggregated by mean"),
             theme = theme_caption),
         "B1_metacells.pdf",
         width = 2.3 * length(all_celltypes) + 1.5, height = 11)

# distribution of |correlations|
n_cor_subsample <- 100000                     # pairs per dataset, plotting only
cor_dist_df <- apply_per_dataset(cor_mats, function(cm, method, ct, smp) {
  v <- abs(cm[upper.tri(cm)])
  data.frame(abs_cor = with_seed(seed, sample(v, min(n_cor_subsample, length(v)))))
})

p_cor_dist <- ggplot(cor_dist_df, aes(abs_cor, colour = condition, group = sample)) +
  geom_density(linewidth = 0.5, adjust = 1.5) +
  facet_grid(method ~ celltype, scales = "free_y",
             labeller = labeller(method = method_labeller)) +
  scale_cond_colour() +
  labs(x = "|Pearson correlation| between gene pairs", y = "Density",
       title = "Correlation distribution across samples, cell types and pooling strategy",
       subtitle = "One curve per sample, coloured by condition",
       caption = fig_caption(
         "Pearson correlation across the metacells within one sample",
         sprintf("%s gene pairs subsampled per dataset for plotting",
                 scales::comma(n_cor_subsample)))) +
  theme(legend.position = "bottom")

save_pdf(p_cor_dist, "B2_correlation_distribution.pdf",
         width = 2.3 * length(all_celltypes) + 1.5, height = 7)

# ---- C GO CO-MEMBERSHIP OF GENE PAIRS ----
go_top_edge_q <- 0.01   # |cor| quantile that defines the "high" edge class

go_annot <- go_shared_matrix(pt_gene_set)
go_edge_df <- apply_per_dataset(cor_mats, function(cm, method, ct, smp)
  go_edge_enrichment(cm, go_annot, go_top_edge_q, n_perm, seed))

go_shift_df <- apply_per_dataset(cor_mats, function(cm, method, ct, smp)
  go_distribution_shift(cm, go_annot, n_perm, seed))

# 2x2 tables, one page per sample, one row per GO term size ceiling
grid_datasets  <- intersect(all_celltypes, c("PT-S1", pooled_celltypes))
grid_max_sizes <- c(50, 100, go_max_term_size)
stat_cols      <- c("a", "b", "c", "d", "n_high", "n_shared", "odds_ratio", "perm_p")

grid_stats <- do.call(rbind, lapply(grid_max_sizes, function(ms) {
  # the widest ceiling is the annotation the rest of the section already ran on
  if (ms == go_max_term_size) {
    dd <- go_edge_df[go_edge_df$method == reference_method & go_edge_df$celltype %in% grid_datasets, ]
    return(data.frame(sample = dd$sample, dataset = as.character(dd$celltype),
                      max_size = ms, dd[, stat_cols], row.names = NULL))
  }
  ann <- go_shared_matrix(pt_gene_set, ms)
  do.call(rbind, lapply(sample_levels, function(smp)
    do.call(rbind, lapply(grid_datasets, function(ds)
      data.frame(sample = smp, dataset = ds, max_size = ms,
                 go_edge_enrichment(cor_mats[[reference_method]][[ds]][[smp]], ann,
                                    go_top_edge_q, n_perm, seed)[, stat_cols],
                 row.names = NULL)))))
}))

as_grid_facets <- function(d) {
  d$dataset <- factor(d$dataset, levels = grid_datasets)
  d$row_lab <- factor(sprintf("GO term size %d-%d", go_min_term_size, d$max_size),
                      levels = sprintf("GO term size %d-%d", go_min_term_size,
                                       grid_max_sizes))
  d
}
grid_stats <- as_grid_facets(grid_stats)
grid_stats$lab <- sprintf("OR = %.2f      permutation p = %s", grid_stats$odds_ratio,
                          format.pval(grid_stats$perm_p, digits = 2,
                                      eps = 1 / (n_perm + 1)))

grid_cells <- as_grid_facets(do.call(rbind, lapply(seq_len(nrow(grid_stats)), function(i) {
  s  <- grid_stats[i, ]
  n  <- as.double(s$a + s$b + s$c + s$d)
  nh <- as.double(s$n_high); ns <- as.double(s$n_shared)
  count    <- c(s$a, s$b, s$c, s$d)
  expected <- c(nh * ns, nh * (n - ns), (n - nh) * ns, (n - nh) * (n - ns)) / n
  data.frame(sample = s$sample, dataset = as.character(s$dataset),
             max_size = s$max_size, x = c(1, 2, 1, 2), y = c(2, 2, 1, 1),
             lab = paste0(scales::comma(count),
                          c(sprintf("\n%.2f\u00d7 expected", count[1] / expected[1]),
                            "", "", "")),
             log2_obs_exp = log2(pmax(count, 0.5) / expected), row.names = NULL)
})))

grid_margins <- as_grid_facets(do.call(rbind, lapply(seq_len(nrow(grid_stats)), function(i) {
  s <- grid_stats[i, ]; n <- as.double(s$a + s$b + s$c + s$d)
  data.frame(sample = s$sample, dataset = as.character(s$dataset),
             max_size = s$max_size, x = c(1, 2, 2.62, 2.62), y = c(0.32, 0.32, 2, 1),
             lab = paste0("n = ", scales::comma(c(s$n_shared, n - s$n_shared,
                                                  s$n_high, n - s$n_high))),
             angle = c(0, 0, 270, 270), row.names = NULL)
})))

grid_axis <- data.frame(
  x = c(1, 2, 0.32, 0.32), y = c(2.64, 2.64, 2, 1),
  lab = c("shared GO term", "no shared term",
          sprintf("top %g%% |cor|", 100 * go_top_edge_q), "remaining pairs"),
  angle = c(0, 0, 90, 90), stringsAsFactors = FALSE)

colour_limit <- max(abs(grid_cells$log2_obs_exp))   # one scale for the whole PDF

grid_page <- function(smp) {
  txt <- function(d, size, ...)
    geom_text(data = d[d$sample == smp, ], aes(x, y, label = lab, angle = angle),
              inherit.aes = FALSE, size = size, ...)
  ggplot(grid_cells[grid_cells$sample == smp, ], aes(x, y)) +
    geom_tile(aes(fill = log2_obs_exp), width = 0.94, height = 0.94) +
    geom_text(aes(label = lab, colour = abs(log2_obs_exp) > 0.55 * colour_limit),
              size = 2.6, lineheight = 0.95, show.legend = FALSE) +
    txt(cbind(grid_axis, sample = smp), 2.4, fontface = "bold", colour = "grey20") +
    txt(grid_margins, 2.1, colour = "grey45") +
    geom_text(data = grid_stats[grid_stats$sample == smp, ],
              aes(1.5, -0.10, label = lab), inherit.aes = FALSE, size = 3,
              fontface = "bold") +
    scale_colour_manual(values = c(`FALSE` = "grey10", `TRUE` = "white")) +
    scale_fill_gradient2(low = "#2166AC", mid = "grey93", high = "#B2182B",
                         midpoint = 0, limits = c(-colour_limit, colour_limit),
                         name = expression(log[2]~"(observed / expected)")) +
    scale_x_continuous(limits = c(0.05, 2.95), expand = c(0, 0)) +
    scale_y_continuous(limits = c(-0.30, 2.95), expand = c(0, 0)) +
    facet_grid(row_lab ~ dataset, switch = "y") +
    theme_void(base_size = 9) +
    theme(legend.position = "right",
          plot.margin   = margin(10, 12, 10, 12),
          panel.border  = element_rect(colour = "grey85", fill = NA, linewidth = 0.3),
          panel.spacing = unit(0.7, "lines"),
          plot.title    = element_text(size = 11, face = "bold", margin = margin(b = 4)),
          plot.subtitle = element_text(size = 8, colour = "grey35", margin = margin(b = 16)),
          plot.caption  = element_text(size = 6.5, colour = "grey45", hjust = 0,
                                       margin = margin(t = 10)),
          strip.text.x      = element_text(size = 8, face = "bold", margin = margin(b = 4)),
          strip.text.y.left = element_text(size = 10, face = "bold", angle = 90,
                                           margin = margin(r = 7))) +
    labs(title = sprintf("GO:BP co-membership of gene pairs within the top %g%% of |cor|, sample %s (%s)",
                         100 * go_top_edge_q, smp,
                         condition_labels[[sample_to_condition[[smp]]]]),
         subtitle = "Shared = two genes are in at least one identical GO:BP term. Colour = log2(observed / expected) under independence. Rows show different GO term size ceiling",
         caption = fig_caption(
           method_labels[[reference_method]],
           "expected counts from the row and column marginals of the same table",
           sprintf("permutation p from %s gene-label permutations",
                   scales::comma(n_perm)),
           "term size filtered on the org.Mm.eg.db scale, ancestor-propagated (GOALL)",
           sprintf("colour scale shared across all %d pages", length(sample_levels))))
}

save_pdf(lapply(sample_levels, grid_page), "C1_GO_2x2_panels.pdf",
         width  = 3.3 * length(grid_datasets) + 1.8,
         height = 3.8 * length(grid_max_sizes) + 1.4)

# odds ratio across the whole edge fraction grid
sweep_q_grid <- seq(0.01, 0.99, by = 0.01)
go_sweep_df  <- apply_per_dataset(cor_mats, function(cm, method, ct, smp)
  go_edge_enrichment_sweep(cm, go_annot, sweep_q_grid))
go_sweep_df  <- go_sweep_df[is.finite(go_sweep_df$odds_ratio), ]

sweep_median <- dplyr::summarise(
  dplyr::group_by(go_sweep_df, method, celltype, top_q),
  odds_ratio = stats::median(odds_ratio), .groups = "drop")

p_sweep <- ggplot(go_sweep_df, aes(top_q, odds_ratio)) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey55") +
  geom_vline(xintercept = go_top_edge_q, linetype = "dotted", colour = "grey35") +
  geom_line(aes(colour = condition, group = sample), alpha = 0.45, linewidth = 0.35) +
  geom_line(data = sweep_median, aes(group = 1), colour = "black", linewidth = 0.8) +
  scale_cond_colour() +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
  facet_grid(method ~ celltype, scales = "free_y",
             labeller = labeller(method = method_labeller)) +
  labs(x = "Fraction of gene pairs classified as high |cor|", y = "Odds ratio",
       title = "Odds ratios displayed for a sweep of the high-|cor| edge fraction",
       subtitle = "Coloured lines = samples, black line = median. Dotted vertical line = the cut-off used in all the other figures",
       caption = fig_caption(
         sprintf("odds ratio (ad)/(bc) of the 2x2 table, threshold grid %g%% to %g%% of pairs",
                 100 * min(sweep_q_grid), 100 * max(sweep_q_grid)),
         sprintf("GO:BP term size %d-%d, ancestor-propagated",
                 go_min_term_size, go_max_term_size))) +
  theme(legend.position = "bottom")

save_pdf(p_sweep, "C2_GO_edge_fraction_sweep.pdf",
         width = 2.3 * length(all_celltypes) + 1.5, height = 7)

# the |cor| densities behind the AUC, one example sample
n_dist_subsample <- 100000       # gene pairs per class and dataset, for plotting only
dist_sample      <- sample_levels[1]

go_dist_long <- do.call(rbind, lapply(all_celltypes, function(ct) {
  pv <- go_pair_vectors(cor_mats[[reference_method]][[ct]][[dist_sample]], go_annot)
  do.call(rbind, lapply(c(TRUE, FALSE), function(is_shared) {
    v <- pv$abs_cor[pv$shared == is_shared]
    data.frame(celltype = ct,
               class = if (is_shared) "shared GO:BP term" else "no shared term",
               abs_cor = with_seed(seed, sample(v, min(n_dist_subsample, length(v)))),
               stringsAsFactors = FALSE)
  }))
}))
go_dist_long$celltype <- factor(go_dist_long$celltype, levels = all_celltypes)

p_dist <- ggplot(go_dist_long, aes(abs_cor, colour = class)) +
  geom_density(linewidth = 0.7, adjust = 1.3) +
  facet_wrap(~ celltype, scales = "free_y") +
  scale_colour_manual(values = c("shared GO:BP term" = "#B2182B",
                                 "no shared term"    = "#4393C3"), name = NULL) +
  labs(x = "|Pearson correlation| between gene pairs", y = "Density",
       title = "Distribution of |cor| for gene pairs with and without a shared GO:BP term",
       subtitle = sprintf("Sample %s (%s), %s", dist_sample,
                          condition_labels[[sample_to_condition[[dist_sample]]]],
                          method_labels[[reference_method]]),
       caption = fig_caption(
         sprintf("GO:BP term size %d-%d, ancestor-propagated; only annotated genes enter the pair set",
                 go_min_term_size, go_max_term_size),
         sprintf("up to %s pairs per class subsampled for plotting",
                 scales::comma(n_dist_subsample)))) +
  theme(legend.position = "bottom")

save_pdf(p_dist, "C3_GO_correlation_distribution.pdf", width = 11, height = 7)

# summary of both co-membership statistics
go_metric_labels <- c(
  odds_ratio = sprintf("Odds ratio\nof co-membership in the top %g%% of |cor|",
                       100 * go_top_edge_q),
  auc = "AUC\nP(|cor| shared pair > |cor| unrelated pair)")

as_metric <- function(d, metric) {
  d$metric <- factor(metric, levels = names(go_metric_labels),
                     labels = go_metric_labels)
  d
}

go_stat_df <- rbind(
  as_metric(data.frame(go_edge_df[, c("method", "celltype", "condition")],
                       value = go_edge_df$odds_ratio), "odds_ratio"),
  as_metric(data.frame(go_shift_df[, c("method", "celltype", "condition")],
                       value = go_shift_df$auc), "auc"))
go_stat_counts <- rbind(
  as_metric(significance_counts(go_edge_df, "odds_ratio", 1), "odds_ratio"),
  as_metric(significance_counts(go_shift_df, "auc", 0.5), "auc"))
go_stat_null <- as_metric(data.frame(null = c(1, 0.5)),
                          names(go_metric_labels))

# rows, not columns: facet_grid frees the y scale per row, and the two statistics
# need different ones while the methods should stay comparable
p_go_stats <- ggplot(go_stat_df, aes(celltype, value, colour = condition)) +
  geom_hline(data = go_stat_null, aes(yintercept = null), linetype = "dashed",
             colour = "grey55") +
  geom_point(position = position_jitter(width = 0.16, height = 0, seed = seed),
             size = 2, alpha = 0.85) +
  geom_text(data = go_stat_counts, aes(celltype, Inf, label = label),
            inherit.aes = FALSE, vjust = 1.5, size = 2.6, colour = "grey35") +
  scale_y_continuous(labels = function(x) {
    d <- 2
    while (d < 4 && anyDuplicated(formatC(x[!is.na(x)], format = "f", digits = d)))
      d <- d + 1
    ifelse(is.na(x), NA_character_, formatC(x, format = "f", digits = d))
  }) +
  scale_cond_colour() +
  facet_grid(metric ~ method, scales = "free_y", switch = "y",
             labeller = labeller(method = method_labeller)) +
  labs(x = NULL, y = NULL,
       title = "GO:BP co-membership of gene pairs, tail test and whole-distribution test",
       subtitle = "One point per sample. Dashed = no association. Numbers above each group: samples with permutation p < 0.05 above the dashed line",
       caption = fig_caption(
         sprintf("odds ratio (ad)/(bc) in the top %g%% of |cor|, GO:BP term size %d-%d, ancestor-propagated",
                 100 * go_top_edge_q, go_min_term_size, go_max_term_size),
         "AUC is the Mann-Whitney statistic over all annotated gene pairs",
         sprintf("%s gene-label permutations per sample and statistic",
                 scales::comma(n_perm)))) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
        legend.position = "bottom",
        strip.placement = "outside",
        strip.text.y.left = element_text(angle = 90))

save_pdf(p_go_stats, "C4_GO_tests_summary.pdf", width = 14, height = 8.5)

# ---- D REPRODUCIBILITY ACROSS SAMPLES ----
rep_top_v <- c(20, 50, 100, 200, 500, 1000)

# on the -z arms every gene has SD ~ 1, so both variants use the -raw ranking
sd_source <- setNames(all_celltypes, all_celltypes)
sd_source[paste0(pool_prefix_mc, "-z")] <- paste0(pool_prefix_mc, "-raw")
sd_source[paste0(pool_prefix_sc, "-z")] <- paste0(pool_prefix_sc, "-raw")

top_sd_genes <- setNames(lapply(metacell_methods, function(method) {
  per_ct <- setNames(lapply(all_celltypes, function(ct)
    lapply(metacells[[method]][[ct]], function(mc) {
      e   <- load_expr_mat(mc, gene_sets[[ct]])
      sds <- setNames(matrixStats::rowSds(e), rownames(e))
      head(names(sort(sds[is.finite(sds) & sds > 1e-8], decreasing = TRUE)),
           max(rep_top_v))
    })), all_celltypes)
  per_ct[all_celltypes] <- per_ct[sd_source[all_celltypes]]
  per_ct
}), metacell_methods)

replicate_pairs <- do.call(rbind, lapply(condition_levels, function(cond) {
  pp <- utils::combn(sample_meta$sample[sample_meta$condition == cond], 2)
  data.frame(condition = cond, sample_a = pp[1, ], sample_b = pp[2, ],
             stringsAsFactors = FALSE)
}))

rep_grid <- expand.grid(top_v = rep_top_v, pair = seq_len(nrow(replicate_pairs)),
                        celltype = all_celltypes, method = metacell_methods,
                        stringsAsFactors = FALSE)

rep_cor_df <- do.call(rbind, lapply(seq_len(nrow(rep_grid)), function(i) {
  g   <- rep_grid[i, ]
  a   <- replicate_pairs$sample_a[g$pair]; b <- replicate_pairs$sample_b[g$pair]
  res <- coexpression_correlation(
    cor_mats[[g$method]][[g$celltype]][[a]], cor_mats[[g$method]][[g$celltype]][[b]],
    intersect(head(top_sd_genes[[g$method]][[g$celltype]][[a]], g$top_v),
              head(top_sd_genes[[g$method]][[g$celltype]][[b]], g$top_v)),
    min_genes = 10L)
  data.frame(method = g$method, celltype = g$celltype,
             condition = replicate_pairs$condition[g$pair],
             pair = paste(a, b, sep = " / "), top_v = g$top_v,
             r = res$r, n_genes = res$n, stringsAsFactors = FALSE)
}))
rep_cor_df$celltype  <- factor(rep_cor_df$celltype, levels = all_celltypes)
rep_cor_df$method    <- factor(rep_cor_df$method, levels = metacell_methods)
rep_cor_df$condition <- factor(rep_cor_df$condition, levels = condition_levels)

p_rep <- ggplot(rep_cor_df, aes(top_v, r, colour = condition, group = pair)) +
  geom_hline(yintercept = 0, linetype = "dotted", colour = "grey55") +
  geom_line(alpha = 0.85, na.rm = TRUE) +
  geom_point(size = 1.6, na.rm = TRUE) +
  ggrepel::geom_text_repel(aes(label = n_genes), size = 2.1, segment.size = 0.2,
                           segment.alpha = 0.5, max.overlaps = Inf, seed = seed,
                           na.rm = TRUE, show.legend = FALSE) +
  scale_x_log10(breaks = rep_top_v) +
  scale_cond_colour() +
  facet_grid(method ~ celltype, labeller = labeller(method = method_labeller)) +
  labs(x = "Number of top-variance genes retained (log scale)",
       y = "Replicate concordance r",
       title = "Pearson cor of co-expression matrices between samples of the same condition",
       subtitle = "One line per replicate pair, colour = condition. Gene sets = intersection of the top-variance genes in both samples",
       caption = fig_caption(
         "r = Pearson correlation between the two upper-triangle correlation vectors",
         "genes ranked by SD across the metacells of the sample, -z arms ranked on their -raw counterpart",
         "point labels = size of the intersected gene set actually used")) +
  theme(legend.position = "bottom",
        axis.text.x = element_text(size = 6, angle = 45, hjust = 1))

save_pdf(p_rep, "D1_replicate_correlation_sweep.pdf",
         width = 2.3 * length(all_celltypes) + 1.5, height = 7)

# all sample pairs, not just replicates
rep_heatmap_top_v <- 500

heat_grid <- expand.grid(sample_b = sample_levels, sample_a = sample_levels,
                         celltype = all_celltypes, stringsAsFactors = FALSE)

heatmap_df <- do.call(rbind, lapply(seq_len(nrow(heat_grid)), function(i) {
  g    <- heat_grid[i, ]
  tops <- top_sd_genes[[reference_method]][[g$celltype]]
  res  <- if (g$sample_a == g$sample_b) list(r = NA_real_, n = NA_real_) else
    coexpression_correlation(cor_mats[[reference_method]][[g$celltype]][[g$sample_a]],
                             cor_mats[[reference_method]][[g$celltype]][[g$sample_b]],
                             intersect(head(tops[[g$sample_a]], rep_heatmap_top_v),
                                       head(tops[[g$sample_b]], rep_heatmap_top_v)))
  data.frame(celltype = g$celltype, sample_a = g$sample_a, sample_b = g$sample_b,
             r = res$r, n_genes = res$n, stringsAsFactors = FALSE)
}))

block      <- as.integer(factor(sample_to_condition[sample_levels],
                                levels = condition_levels))
axis_pos   <- setNames(seq_along(sample_levels) + 0.15 * (block - 1), sample_levels)
axis_label <- setNames(paste0(sample_levels, "  |  ", sample_to_batch[sample_levels]),
                       sample_levels)
axis_col   <- unname(cond_colours[sample_to_condition[sample_levels]])

heatmap_df$x   <- axis_pos[heatmap_df$sample_a]
heatmap_df$y   <- axis_pos[heatmap_df$sample_b]
heatmap_df$rel <- stats::ave(heatmap_df$r, heatmap_df$celltype, FUN = function(v) {
  rng <- diff(range(v, na.rm = TRUE))
  if (!is.finite(rng) || rng == 0) 0.5 else (v - min(v, na.rm = TRUE)) / rng
})
# blank facet levels pad the first row so that base and pooled cell types line up
pad_levels <- strrep(" ", seq_len(max(0L, length(pooled_celltypes) -
                                        length(base_celltypes))))
heatmap_df$celltype <- factor(heatmap_df$celltype, levels = all_celltypes)

# base and pooled cell types as two rows; the shorter row leaves its last cell empty
heat_row <- function(d) ggplot(d, aes(x, y, fill = rel)) +
  geom_tile(width = 1, height = 1, colour = "white", linewidth = 0.25) +
  # "+ 0" turns a rounded -0 into 0, otherwise sprintf prints "-0"
  geom_text(aes(label = ifelse(is.na(r), "", sprintf("%.0f", round(100 * r) + 0)),
                colour = !is.na(rel) & rel > 0.6), size = 2.3, show.legend = FALSE) +
  scale_fill_gradient(low = "#F7FBFF", high = "#08306B", na.value = "grey92",
                      breaks = c(0, 1), labels = c("lowest", "highest"),
                      name = "within panel") +
  scale_colour_manual(values = c(`FALSE` = "grey15", `TRUE` = "white"), guide = "none") +
  scale_x_continuous(breaks = axis_pos, labels = axis_label,
                     expand = expansion(add = 0.6)) +
  scale_y_reverse(breaks = axis_pos, labels = axis_label,
                  expand = expansion(add = 0.6)) +
  facet_wrap(~ celltype, ncol = max(length(base_celltypes), length(pooled_celltypes))) +
  coord_fixed() +
  labs(x = NULL, y = NULL) +
  theme_qc(base_size = 9) +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 5.5,
                                   colour = axis_col),
        axis.text.y = element_text(size = 5.5, colour = axis_col),
        axis.ticks  = element_blank())

save_pdf(patchwork::wrap_plots(
  heat_row(droplevels(heatmap_df[heatmap_df$celltype %in% base_celltypes, ])),
  heat_row(droplevels(heatmap_df[heatmap_df$celltype %in% pooled_celltypes, ])),
  ncol = 1, guides = "collect") +
    patchwork::plot_annotation(
      title = "Pearson cor of co-expression matrices between all sample pairs",
      subtitle = sprintf("Panel numbers = Pearson cor x 100, %s. Samples are grouped by condition",
                         method_labels[[reference_method]]),
      caption = fig_caption(
        sprintf("gene set = top %d by SD in A intersected with top %d in B (%d-%d genes per cell)",
                rep_heatmap_top_v, rep_heatmap_top_v,
                min(heatmap_df$n_genes, na.rm = TRUE),
                max(heatmap_df$n_genes, na.rm = TRUE)),
        "colour is scaled within each panel",
        "CR is ~75% confounded with batch TX125"),
      theme = theme_caption),
  "D2_sample_correlation_heatmap.pdf", width = 14, height = 8)

# ---- E SAMPLE-LEVEL PCA ----
plot_sample_pca <- function(pca, title, k_sep = 2L) {
  var_exp <- pca$sdev^2 / sum(pca$sdev^2)
  meta    <- sample_meta[match(rownames(pca$x), sample_meta$sample), ]
  coords  <- pca$x[, seq_len(min(k_sep, ncol(pca$x))), drop = FALSE]
  df      <- merge(data.frame(sample = rownames(pca$x), pca$x[, 1:2, drop = FALSE],
                              stringsAsFactors = FALSE), sample_meta, by = "sample")
  df$condition <- factor(df$condition, levels = condition_levels)
  
  ggplot(df, aes(PC1, PC2, colour = condition, shape = batch)) +
    geom_point(size = 3) +
    ggrepel::geom_text_repel(aes(label = sample), size = 2.8, colour = "grey25",
                             min.segment.length = 0.2, seed = seed,
                             show.legend = FALSE) +
    scale_cond_colour() +
    scale_shape_manual(values = batch_shapes, name = "Batch") +
    labs(title = title,
         subtitle = sprintf("separation in PC1-%d (mean between / mean within distance): condition %.2f, batch %.2f",
                            ncol(coords), separation_ratio(coords, meta$condition),
                            separation_ratio(coords, meta$batch)),
         x = sprintf("PC1 (%.1f%% of variance)", 100 * var_exp[1]),
         y = sprintf("PC2 (%.1f%% of variance)", 100 * var_exp[2]))
}

pca_pages <- lapply(all_celltypes, function(ct) {
  panels <- list(
    expression = plot_sample_pca(
      sample_pca(metacells[[reference_method]][[ct]], "expression", gene_sets[[ct]]),
      "mean expression space"),
    structure = plot_sample_pca(
      sample_pca(cor_mats[[reference_method]][[ct]], "structure"),
      "co-expression structure space"))
  patchwork::wrap_plots(panels, ncol = 2, guides = "collect") &
    theme(legend.position = "bottom", plot.margin = margin(4, 6, 4, 6))
})
pca_pages <- lapply(seq_along(pca_pages), function(i)
  pca_pages[[i]] + patchwork::plot_annotation(
    title = sprintf("%s  |  %s", all_celltypes[i], method_labels[[reference_method]]),
    caption = fig_caption(
      "expression space: mean of the data layer per gene, features scaled before the PCA",
      "structure space: upper triangle of the correlation matrix, unscaled"),
    theme = theme_caption + theme(plot.title = element_text(hjust = 0.5))))

save_pdf(pca_pages, "E1_sample_PCA.pdf", width = 12, height = 6)

# ---- F PER-GENE CONNECTIVITY MODEL ----
conn_p        <- 0.01
conn_scalings <- list(
  absolute = identity,                                   # connectivity as computed
  relative = function(x) sweep(x, 2, colMeans(x), "/"))  # divided by the sample mean

conn_celltypes    <- paste0(c(pool_prefix_sc, pool_prefix_mc), "-raw")
conn_term_levels  <- c("age", "diet", "interaction")
conn_term_labels  <- c(age         = "Age (old - young, averaged over diet)",
                       diet        = "Diet (CR - AL, averaged over age)",
                       interaction = "Age x Diet interaction")
conn_scale_labels <- c(absolute = "absolute: mean |atanh(r)|",
                       relative = "relative: mean |atanh(r)| divided by the sample mean")

# +-0.5 effect coding; trend and robust both matter at 4 residual df
conn_fit <- function(conn) {
  cond <- sample_to_condition[colnames(conn)]
  des  <- stats::model.matrix(~ age_c * diet_c, data = data.frame(
    age_c  = ifelse(sub(".*_", "", cond) == "old", 0.5, -0.5),
    diet_c = ifelse(sub("_.*", "", cond) == "CR",  0.5, -0.5)))
  dimnames(des) <- list(colnames(conn), c("intercept", "age", "diet", "interaction"))
  fit <- limma::eBayes(limma::lmFit(conn, des), trend = TRUE, robust = TRUE)
  do.call(rbind, lapply(conn_term_levels, function(cf) {
    tt <- limma::topTable(fit, coef = cf, number = Inf, sort.by = "none")
    data.frame(gene = rownames(tt), term = cf, slope = tt$logFC, t = tt$t,
               p = tt$P.Value, row.names = NULL, stringsAsFactors = FALSE)
  }))
}

conn_n_hvg    <- 2000
conn_gene_set <- hvg_sets[[as.character(conn_n_hvg)]]$final
conn_raw <- setNames(lapply(conn_celltypes, function(ct)
  connectivity_matrix(lapply(cor_mats[[reference_method]][[ct]][sample_levels],
                             function(cm) {
                               g <- intersect(conn_gene_set, rownames(cm))
                               cm[g, g, drop = FALSE]
                             }))),
  conn_celltypes)

conn_stats <- do.call(rbind, lapply(conn_celltypes, function(ct)
  do.call(rbind, lapply(names(conn_scalings), function(sc)
    data.frame(celltype = ct, scaling = sc,
               conn_fit(conn_scalings[[sc]](conn_raw[[ct]])),
               row.names = NULL, stringsAsFactors = FALSE)))))

conn_stats$celltype <- factor(conn_stats$celltype, levels = conn_celltypes)
conn_stats$scaling  <- factor(conn_stats$scaling, levels = names(conn_scalings),
                              labels = conn_scale_labels[names(conn_scalings)])
conn_stats$term     <- factor(conn_stats$term, levels = conn_term_levels,
                              labels = conn_term_labels[conn_term_levels])
arm_levels <- as.vector(outer(sub("^(.)", "\\U\\1", names(conn_scalings), perl = TRUE),
                              conn_celltypes, function(sc, ct) sprintf("%s; %s", ct, sc)))
conn_stats$arm <- factor(
  sprintf("%s; %s", conn_stats$celltype,
          sub("^(.)", "\\U\\1", sub(":.*", "", conn_stats$scaling), perl = TRUE)),
  levels = arm_levels)

conn_counts <- dplyr::summarise(
  dplyr::group_by(conn_stats, celltype, scaling, arm, term),
  n_genes = dplyr::n(), n_sig = sum(p < conn_p), exp_sig = conn_p * dplyr::n(),
  label = sprintf("%d / %s at p < %g\n%.0f expected by chance",
                  sum(p < conn_p), scales::comma(dplyr::n()), conn_p,
                  conn_p * dplyr::n()),
  .groups = "drop")

# per-sample metrics behind the model, one facet each
mol_panel_celltype <- conn_celltypes[1]

# mean expression and variance the correlations were computed from
conn_var_df <- do.call(rbind, lapply(sample_levels, function(s) {
  e <- load_expr_mat(metacells[[reference_method]][[mol_panel_celltype]][[s]],
                     rownames(conn_raw[[mol_panel_celltype]]))
  data.frame(sample = s, variance = matrixStats::rowVars(e), expr = rowMeans(e),
             row.names = NULL, stringsAsFactors = FALSE)
}))
conn_var_df <- conn_var_df[is.finite(conn_var_df$variance) & conn_var_df$variance > 0, ]

lib_df <- do.call(rbind, lapply(c("pre", "post"), function(ds) {
  d <- if (ds == "pre") lib_prenorm_df else lib_all_df
  rbind(data.frame(key = paste0("lib_", ds, "_all"), d[, c("sample", "lib")]),
        data.frame(key = paste0("lib_", ds, "_pt"),
                   d[d$cell %in% unname(cells_kept), c("sample", "lib")]))
}))
lib_df <- lib_df[lib_df$lib > 0 & lib_df$sample %in% sample_levels, ]

mol_metric_levels <- c(
  lib_pre_all  = sprintf("Library size (log10)\nbefore depth norm.\nall %d cell types",
                         length(unique(lib_all_df$celltype))),
  lib_pre_pt   = "Library size (log10)\nbefore depth norm.\nPT only, equalised",
  lib_post_all = sprintf("Library size (log10)\nafter depth norm.\nall %d cell types",
                         length(unique(lib_all_df$celltype))),
  lib_post_pt  = "Library size (log10)\nafter depth norm.\nPT only, equalised",
  singleton    = "Molecules with\nexactly one read",
  expr         = "Mean expression\nper gene",
  variance     = "Variance \nper gene",
  conn         = "Connectivity (mean |atanh(r)|)\nper gene")

as_metric_df <- function(sample, metric, value)
  data.frame(sample = as.character(sample), metric = metric, value = value,
             row.names = NULL, stringsAsFactors = FALSE)

mol_dist <- rbind(
  as_metric_df(lib_df$sample, lib_df$key, log10(lib_df$lib)),
  as_metric_df(conn_var_df$sample, "expr",     conn_var_df$expr),
  as_metric_df(conn_var_df$sample, "variance", conn_var_df$variance),
  as_metric_df(rep(colnames(conn_raw[[mol_panel_celltype]]),
                   each = nrow(conn_raw[[mol_panel_celltype]])), "conn",
               as.numeric(conn_raw[[mol_panel_celltype]])))

add_metric_meta <- function(d) {
  d$metric    <- factor(d$metric, levels = names(mol_metric_levels),
                        labels = mol_metric_levels)
  d$sample    <- factor(d$sample, levels = sample_levels)
  d$condition <- factor(sample_to_condition[as.character(d$sample)],
                        levels = condition_levels)
  d
}
mol_dist  <- add_metric_meta(mol_dist)
mol_point <- add_metric_meta(as_metric_df(mol_stats$sample, "singleton",
                                          mol_stats$singleton))

mol_log_metrics <- mol_metric_levels[c("lib_pre_all", "lib_pre_pt", "lib_post_all",
                                       "lib_post_pt", "variance")]
mol_med <- dplyr::summarise(
  dplyr::group_by(mol_dist, metric, sample),
  med  = stats::median(value),
  mean = if (metric[1] %in% mol_log_metrics) log10(mean(10^value)) else mean(value),
  .groups = "drop")

mol_reads_lab <- do.call(rbind, lapply(
  c("lib_pre_all", "lib_pre_pt", "lib_post_all", "lib_post_pt"), function(k)
    data.frame(metric = mol_metric_levels[[k]],
               sample = as.character(mol_stats$sample),
               lab = sprintf("%.0fk", (if (grepl("_pre_", k)) mol_stats$reads_pc
                                       else mol_target) / 1000),
               row.names = NULL, stringsAsFactors = FALSE)))
mol_reads_lab$metric <- factor(mol_reads_lab$metric, levels = mol_metric_levels)
mol_reads_lab$sample <- factor(mol_reads_lab$sample, levels = sample_levels)

p_mol_overview <- ggplot(mol_dist, aes(sample, value)) +
  geom_violin(aes(fill = condition), scale = "width", linewidth = 0.3,
              colour = "grey45") +
  geom_errorbar(data = mol_med, aes(x = sample, ymin = med, ymax = med),
                inherit.aes = FALSE, width = 0.8, linewidth = 0.4,
                linetype = "dashed", colour = "grey15") +
  geom_errorbar(data = mol_med, aes(x = sample, ymin = mean, ymax = mean),
                inherit.aes = FALSE, width = 0.8, linewidth = 0.4, colour = "grey15") +
  geom_point(data = mol_point, aes(colour = condition), shape = 16, size = 3.2,
             show.legend = FALSE) +
  geom_text(data = mol_reads_lab, aes(sample, Inf, label = lab), inherit.aes = FALSE,
            vjust = 1.4, size = 1.9, colour = "grey35") +
  facet_wrap(~ metric, nrow = 1, scales = "free_y") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.12))) +
  scale_x_discrete(labels = function(s) paste0(s, "  |  ", sample_to_batch[s])) +
  scale_cond_fill() +
  scale_cond_colour() +
  labs(x = NULL, y = NULL,
       title = "Input for the linear connectivity model",
       subtitle = sprintf("Library size plots showing reads per cell above each violin; Expression, variance and connectivity per gene from %s; solid line = mean, dashed line = median",
                          mol_panel_celltype),
       caption = fig_caption(
         method_labels[[reference_method]],
         "library size from the counts layer, pre-normalisation values from the object before the CellRanger aggr depth normalisation, cells matched by barcode",
         sprintf("aggr target %s reads per cell",
                 scales::comma(round(mol_target))),
         "singleton rate from molecule_info.h5",
         "expression and variance from data layer, library size divided out per metacell before the log",
         "connectivity = mean |atanh(Pearson r)| over all other genes")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6),
        legend.position = "bottom")

save_pdf(p_mol_overview, "F1_model_input.pdf",
         width = 2.1 * length(mol_metric_levels) + 2, height = 6)

# coefficients
conn_lab <- dplyr::slice_min(dplyr::group_by(conn_stats, celltype, scaling, term),
                             order_by = p, n = 5, with_ties = FALSE)

p_conn_volcano <- ggplot(conn_stats, aes(slope, -log10(p))) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey55") +
  geom_hline(yintercept = -log10(conn_p), linetype = "dotted", colour = "#B2182B",
             linewidth = 0.4) +
  geom_point(aes(colour = p < conn_p), size = 0.6, alpha = 0.5) +
  ggrepel::geom_text_repel(data = conn_lab, aes(label = gene), size = 2,
                           max.overlaps = Inf, seed = seed, colour = "grey20") +
  geom_text(data = conn_counts, aes(-Inf, Inf, label = label), inherit.aes = FALSE,
            hjust = -0.03, vjust = 1.15, size = 2.2, colour = "grey35",
            lineheight = 0.85) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.22))) +
  scale_colour_manual(values = c(`FALSE` = "grey70", `TRUE` = "#B2182B"),
                      labels = c(sprintf("p >= %g", conn_p),
                                 sprintf("p < %g", conn_p)), name = NULL) +
  facet_grid(term ~ arm, scales = "free",
             labeller = labeller(
               arm  = function(x) paste0(x, "\n", ifelse(grepl("Absolute$", x),
                                                         "change in mean |atanh(r)|",
                                                         "change in multiples of the sample mean")),
               term = label_wrap_gen(22))) +
  labs(x = "Coefficient", y = expression(-log[10]~"p"),
       title = "Linear connectivity model coefficients",
       subtitle = sprintf("Absolute = raw connectivity values used for lm; relative = connectivities divided by the sample mean across all genes; dotted line = p < %g",
                          conn_p),
       caption = fig_caption(
         "~ age * diet with +-0.5 effect coding, limma lmFit without batch term, eBayes(trend = TRUE, robust = TRUE), 8 samples and 4 residual df before moderation",
         "SE(interaction) is twice SE(age) and SE(diet)",
         "labels = the 5 smallest p per panel")) +
  theme(legend.position = "bottom", axis.text.x = element_text(size = 6),
        strip.text.x = element_text(size = 8, margin = margin(5, 3, 5, 3)))

save_pdf(p_conn_volcano, "F2_connectivity_coefficients.pdf",
         width = 3.2 * length(conn_celltypes) * length(conn_scalings) + 2, height = 9)

# GSEA on the ranked connectivity coefficients
gsea_min_size   <- 8        # on universe level
gsea_show_n     <- 5        # terms per cluster, picked by enrichplot on the nominal p
gsea_term_short <- c(age = "Age", diet = "Diet", interaction = "Age x Diet")
gsea_rank_short <- c(t = "moderated t", slope = "raw coef.")

run_gsea <- function(stat) {
  stat <- sort(stat[is.finite(stat)], decreasing = TRUE)
  with_seed(seed, suppressWarnings(clusterProfiler::GSEA(
    geneList = stat, TERM2GENE = go_term2gene, TERM2NAME = go_term2name,
    minGSSize = gsea_min_size, maxGSSize = length(stat) - gsea_min_size,
    exponent = 1, eps = 0, seed = TRUE, pvalueCutoff = 1, pAdjustMethod = "BH",
    verbose = FALSE)))
}

gsea_objs <- list()
for (ct in conn_celltypes) for (sc in names(conn_scalings)) {
  d <- conn_stats[conn_stats$celltype == ct &
                    conn_stats$scaling == conn_scale_labels[[sc]], ]
  for (cf in conn_term_levels) {
    dd <- d[d$term == conn_term_labels[[cf]], ]
    for (rs in names(gsea_rank_short))
      gsea_objs[[paste(ct, sc, cf, rs, sep = "|")]] <-
        run_gsea(setNames(dd[[rs]], dd$gene))
  }
}

gsea_nes_lim  <- max(abs(unlist(lapply(gsea_objs, function(o) o@result$NES))),
                     na.rm = TRUE)
gsea_padj_lim <- max(-log10(pmax(unlist(lapply(gsea_objs, function(o)
  o@result$p.adjust)), 1e-10)), na.rm = TRUE)

cc_size_range <- c(1.5, 5.5)   # point size, mm
cc_row_mm     <- 10            # row height on the page with the most terms
cc_nudge      <- (max(cc_size_range) / 2 + 1.8) / cc_row_mm
cc_label_n    <- 55            # term name truncation

cc_cols <- expand.grid(rank = names(gsea_rank_short), term = conn_term_levels,
                       stringsAsFactors = FALSE)          # ranking varies fastest
cc_col_levels <- sprintf("%s\n%s", gsea_term_short[cc_cols$term],
                         gsea_rank_short[cc_cols$rank])

gsea_df <- do.call(rbind, lapply(names(gsea_objs), function(k) {
  r <- gsea_objs[[k]]@result
  if (!nrow(r)) return(NULL)
  key <- strsplit(k, "|", fixed = TRUE)[[1]]
  data.frame(ct = key[1], sc = key[2], term = key[3], rank = key[4],
             ID = r$ID, Description = r$Description, NES = r$NES,
             p = r$pvalue, p_adj = r$p.adjust,
             n_core = lengths(strsplit(r$core_enrichment, "/", fixed = TRUE)),
             row.names = NULL, stringsAsFactors = FALSE)
}))
gsea_df$col <- factor(sprintf("%s\n%s", gsea_term_short[gsea_df$term],
                              gsea_rank_short[gsea_df$rank]),
                      levels = cc_col_levels)

# one page per pooling arm, x axis = contrast x ranking
cc_page <- function(ct, sc) {
  d   <- gsea_df[gsea_df$ct == ct & gsea_df$sc == sc, ]
  sel <- unique(unlist(lapply(split(d, d$col), function(x)
    head(x$ID[order(x$p)], gsea_show_n))))
  d   <- d[d$ID %in% sel, ]
  
  ord  <- names(sort(tapply(d$p_adj, d$ID, min), decreasing = TRUE))  # best at top
  term <- d$Description[match(ord, d$ID)]
  term <- make.unique(ifelse(nchar(term) > cc_label_n,
                             paste0(substr(term, 1, cc_label_n - 1), "\u2026"), term))
  d$ID <- factor(d$ID, levels = ord, labels = term)
  
  ggplot(d, aes(col, ID)) +
    geom_point(aes(size = -log10(pmax(p_adj, 1e-10)), fill = NES),
               shape = 21, colour = "grey30", stroke = 0.3) +
    geom_text(aes(label = n_core), nudge_y = cc_nudge, size = 1.9,
              colour = "grey25") +
    scale_x_discrete(drop = FALSE) +
    scale_y_discrete(expand = expansion(add = c(0.6, 0.9))) +
    scale_fill_gradient2(low = "#2166AC", mid = "grey93", high = "#B2182B",
                         midpoint = 0, limits = c(-gsea_nes_lim, gsea_nes_lim),
                         name = "NES") +
    scale_size_continuous(range = cc_size_range, limits = c(0, gsea_padj_lim),
                          name = expression(-log[10]~"adj. p")) +
    labs(x = NULL, y = NULL,
         title = sprintf("%s; %s", ct, sub("^(.)", "\\U\\1", sc, perl = TRUE)),
         subtitle = sprintf("Top %d GO:BP terms per contrast x ranking, selected on the nominal p; rows ordered by the smallest adjusted p across columns",
                            gsea_show_n),
         caption = fig_caption(
           sprintf("clusterProfiler::GSEA per contrast and ranking, GO:BP term size %d-%d (ancestor-propagated), additional universe filter min_universe_genes = %d, BH adjustment",
                   go_min_term_size, go_max_term_size, gsea_min_size),
           "parent-child GO relations are not collapsed",
           "t = slope / (moderated SD); numbers above the dots = genes in the leading edge")) +
    theme(axis.text.y   = element_text(size = 6.5),
          axis.text.x   = element_text(size = 7, lineheight = 1.1),
          plot.caption  = element_text(size = 6.5))
}

cc_pages   <- unlist(lapply(conn_celltypes, function(ct)
  lapply(names(conn_scalings), function(sc) cc_page(ct, sc))), recursive = FALSE)
cc_n_terms <- max(vapply(cc_pages, function(p) nlevels(p$data$ID), integer(1)))

save_pdf(cc_pages, "F3_GSEA.pdf", width = 10, height = cc_n_terms * cc_row_mm / 25.4 + 3.5)

# ---- G TERM-LEVEL MODULE CONNECTIVITY ----
blk_min_genes   <- 10    # genes of a term inside the co-expression universe
blk_n_bg        <- 50000 # decoy sets per size, drawn once and reused per sample
blk_top_n       <- 10    # terms per page
blk_max_overlap <- 0.7   # |A n B| / min(|A|,|B|), catches propagated parent/child
blk_max_jaccard <- 0.5   # |A n B| / |A u B|, catches overlap without nesting
blk_rep_n_label <- 10    # labelled terms

blk_celltype <- paste0(pool_prefix_sc, "-raw")
blk_universe <- Reduce(intersect, lapply(
  cor_mats[[reference_method]][[blk_celltype]][sample_levels], rownames))
blk_hit  <- go_term2gene$SYMBOL %in% blk_universe
blk_sets <- lapply(split(go_term2gene$SYMBOL[blk_hit], go_term2gene$GO[blk_hit]),
                   unique)
blk_sets <- blk_sets[lengths(blk_sets) >= blk_min_genes]

# atanh|r| once per sample: the transform does not depend on the set and the
# background alone needs several hundred thousand block means
blk_A <- lapply(cor_mats[[reference_method]][[blk_celltype]][sample_levels],
                function(cm) {
                  a <- atanh(pmin(abs(cm), 0.999))
                  diag(a) <- 0
                  a
                })

# mean atanh|r| over the m(m-1)/2 within-block pairs
blk_mean <- function(A, g) sum(A[g, g]) / (as.double(length(g)) * (length(g) - 1))

# decoys drawn once per size and reused across samples
blk_sizes  <- sort(unique(lengths(blk_sets)))
blk_decoys <- setNames(lapply(blk_sizes, function(m)
  with_seed(seed + m, lapply(seq_len(blk_n_bg), function(i)
    sample(blk_universe, m)))), as.character(blk_sizes))

blk_bg_val <- setNames(lapply(blk_sizes, function(m)
  vapply(sample_levels, function(s)
    vapply(blk_decoys[[as.character(m)]], blk_mean, numeric(1), A = blk_A[[s]]),
    numeric(blk_n_bg))), as.character(blk_sizes))

blk_bg <- do.call(rbind, lapply(sample_levels, function(s)
  do.call(rbind, lapply(blk_sizes, function(m) {
    v <- blk_bg_val[[as.character(m)]][, s]
    data.frame(sample = s, n_genes = m, bg_mean = mean(v), bg_sd = stats::sd(v),
               row.names = NULL, stringsAsFactors = FALSE)
  }))))

blk_df <- do.call(rbind, lapply(names(blk_sets), function(go)
  data.frame(GO = go, n_genes = length(blk_sets[[go]]), sample = sample_levels,
             obs = vapply(blk_A, blk_mean, numeric(1), g = blk_sets[[go]]),
             row.names = NULL, stringsAsFactors = FALSE)))
blk_df <- merge(blk_df, blk_bg, by = c("sample", "n_genes"))
blk_df$delta <- blk_df$obs - blk_df$bg_mean
blk_df$z     <- blk_df$delta / blk_df$bg_sd
blk_df$condition <- factor(sample_to_condition[blk_df$sample], levels = condition_levels)
blk_df$age   <- factor(sub(".*_", "", as.character(blk_df$condition)),
                       levels = c("young", "old"))
blk_df$diet  <- factor(sub("_.*", "", as.character(blk_df$condition)),
                       levels = c("AL", "CR"))
blk_df$sample <- factor(blk_df$sample, levels = sample_levels)

blk_cond <- dplyr::summarise(
  dplyr::group_by(blk_df, GO, n_genes, condition, age, diet),
  obs = mean(obs), bg_mean = mean(bg_mean), z = mean(z), .groups = "drop")

blk_wide <- tapply(blk_cond$z, list(blk_cond$GO, as.character(blk_cond$condition)),
                   function(x) x[1])

# pairs drive the segments drawn on the G1 pages
blk_contrasts <- list(
  age_AL      = list(pairs = list(c("AL_young", "AL_old")),
                     title = "Age effect within AL"),
  age_CR      = list(pairs = list(c("CR_young", "CR_old")),
                     title = "Age effect within CR"),
  diet_young  = list(pairs = list(c("AL_young", "CR_young")),
                     title = "Diet effect within young"),
  diet_old    = list(pairs = list(c("AL_old", "CR_old")),
                     title = "Diet effect within old"),
  interaction = list(pairs = list(c("AL_young", "AL_old"), c("CR_young", "CR_old")),
                     w = c(AL_young = 1, AL_old = -1, CR_young = -1, CR_old = 1),
                     title = "Age x Diet interaction: \u0394z(CR) - \u0394z(AL)"))

blk_weights <- function(cc) {
  if (!is.null(cc$w)) return(cc$w[condition_levels])
  setNames(as.numeric(condition_levels == cc$pairs[[1]][2]) -
             as.numeric(condition_levels == cc$pairs[[1]][1]), condition_levels)
}

blk_rank <- do.call(rbind, lapply(names(blk_contrasts), function(k)
  data.frame(contrast = k, GO = rownames(blk_wide),
             change = as.numeric(blk_wide[, condition_levels] %*% blk_weights(blk_contrasts[[k]])),
             row.names = NULL, stringsAsFactors = FALSE)))

blk_cond_samples <- setNames(lapply(condition_levels, function(cd)
  sample_levels[sample_to_condition[sample_levels] == cd]), condition_levels)

# every decoy set run through the same condition-mean difference as the real terms -> p value
blk_null_delta <- setNames(lapply(blk_sizes, function(m) {
  v  <- blk_bg_val[[as.character(m)]]
  z  <- sweep(sweep(v, 2, colMeans(v), "-"), 2, matrixStats::colSds(v), "/")
  zc <- vapply(condition_levels, function(cd)
    rowMeans(z[, blk_cond_samples[[cd]], drop = FALSE]), numeric(blk_n_bg))
  vapply(blk_contrasts, function(cc) as.numeric(zc %*% blk_weights(cc)),
         numeric(blk_n_bg))
}), as.character(blk_sizes))

blk_rank$p <- mapply(function(go, k, ch) {
  nd <- blk_null_delta[[as.character(length(blk_sets[[go]]))]][, k]
  (sum(abs(nd) >= abs(ch)) + 1) / (blk_n_bg + 1)
}, blk_rank$GO, blk_rank$contrast, blk_rank$change)
blk_rank$p_adj <- stats::ave(blk_rank$p, blk_rank$contrast,
                             FUN = function(p) stats::p.adjust(p, "BH"))

# a propagated child sits fully inside its parent while Jaccard stays low, so the
# overlap coefficient is the criterion that catches nesting
blk_redundant <- function(a, b) {
  n <- length(intersect(a, b))
  n / min(length(a), length(b)) > blk_max_overlap ||
    n / length(union(a, b))     > blk_max_jaccard
}

# greedy walk down the ranked list, keeping a term only if it is distinct from
# every better ranked term already kept
blk_pick_top <- function(go_ranked, n) {
  keep <- character(0)
  for (go in go_ranked) {
    if (length(keep) >= n) break
    if (!any(vapply(keep, function(kk)
      blk_redundant(blk_sets[[go]], blk_sets[[kk]]), logical(1))))
      keep <- c(keep, go)
  }
  keep
}

blk_line_colours <- c(AL = unname(cond_colours[["AL_young"]]),
                      CR = unname(cond_colours[["CR_young"]]))

blk_page <- function(k) {
  cc  <- blk_contrasts[[k]]
  r   <- blk_rank[blk_rank$contrast == k, ]
  r   <- r[order(-abs(r$change)), ]
  sel <- r[match(blk_pick_top(r$GO, blk_top_n), r$GO), ]
  
  term <- go_term2name$TERM[match(sel$GO, go_term2name$GO)]
  term <- ifelse(nchar(term) > 42, paste0(substr(term, 1, 41), "..."), term)
  lab  <- setNames(sprintf("%s\nm = %d  |  \u0394z = %+.2f", make.unique(term),
                           lengths(blk_sets)[sel$GO], sel$change), sel$GO)
  
  as_page <- function(d) {
    d$GO <- factor(d$GO, levels = sel$GO, labels = lab[sel$GO])
    d
  }
  d <- as_page(blk_df[blk_df$GO %in% sel$GO, ])
  m <- as_page(blk_cond[blk_cond$GO %in% sel$GO, ])
  
  # the ranked edge: a diet line for the age contrasts, the vertical gap between
  # the two diet lines at one age for the diet contrasts
  pick <- function(cond) m$obs[match(paste(levels(d$GO), cond),
                                     paste(as.character(m$GO), m$condition))]
  seg <- do.call(rbind, lapply(cc$pairs, function(pr) data.frame(
    GO   = factor(levels(d$GO), levels = levels(d$GO)),
    x    = as.numeric(factor(sub(".*_", "", pr[1]), levels = c("young", "old"))),
    xend = as.numeric(factor(sub(".*_", "", pr[2]), levels = c("young", "old"))),
    y    = pick(pr[1]), yend = pick(pr[2]))))
  
  seg_lab <- paste(vapply(cc$pairs, function(pr)
    sprintf("%s -> %s", condition_labels[[pr[1]]], condition_labels[[pr[2]]]),
    character(1)), collapse = ", ")
  
  ggplot(d, aes(age, obs)) +
    geom_line(data = m, aes(y = bg_mean, group = diet), colour = "grey70",
              linetype = "22", linewidth = 0.4) +
    geom_segment(data = seg, aes(x = x, xend = xend, y = y, yend = yend),
                 inherit.aes = FALSE, colour = "grey15", linewidth = 1.5,
                 alpha = 0.55) +
    geom_line(data = m, aes(group = diet, colour = diet), linewidth = 0.5) +
    geom_point(aes(fill = condition), shape = 21, colour = "grey30", size = 2,
               stroke = 0.3,
               position = position_jitter(width = 0.07, height = 0, seed = seed)) +
    facet_wrap(~ GO, scales = "free_y", nrow = 2) +
    scale_colour_manual(values = blk_line_colours, name = "Diet (condition mean)") +
    scale_cond_fill() +
    labs(x = NULL, y = "Block connectivity (mean |atanh(r)|)",
         title = sprintf("%s: top %d GO:BP terms by |\u0394z| in background-corrected block connectivity",
                         cc$title, blk_top_n),
         subtitle = sprintf("Points = the 8 samples, coloured lines = condition means, dashed grey = size-matched random background of the same sample. %s. Panels ordered by |\u0394z|",
                            if (length(cc$pairs) > 1)
                              sprintf("Bold segments = the two age slopes whose difference is ranked (%s)", seg_lab)
                            else
                              sprintf("Bold segment = the ranked contrast (%s)", seg_lab)),
         caption = fig_caption(
           blk_celltype, method_labels[[reference_method]],
           sprintf("%s candidate terms with >= %d of the %s universe genes, GO:BP size %d-%d on the org.Mm.eg.db scale, ancestor-propagated",
                   scales::comma(length(blk_sets)), blk_min_genes,
                   scales::comma(length(blk_universe)), go_min_term_size,
                   go_max_term_size),
           "term connectivity = mean atanh|r| over the within-term gene pairs",
           paste0(
             sprintf("background = %s random sets per size, drawn once and reused across samples",
                     scales::comma(blk_n_bg)),
             "\nz = (observed - background mean) / background SD of the same sample and set size, used for the ranking"),
           sprintf("redundancy filter on the ranked list: overlap coefficient > %.2f or Jaccard > %.2f drops a term",
                   blk_max_overlap, blk_max_jaccard))) +
    theme(legend.position = "bottom",
          strip.text = element_text(size = 6.5, lineheight = 1.2))
}

save_pdf(lapply(names(blk_contrasts), blk_page), "G1_block_connectivity.pdf",
         width = 15, height = 8, encoding = "Greek.enc")

# connectivity dotplot
gdot_top_n      <- 5
gdot_order      <- c("diet_young", "diet_old", "age_AL", "age_CR", "interaction")
gdot_col_labels <- c(diet_young  = "Diet\nwithin young", diet_old = "Diet\nwithin old",
                     age_AL      = "Age\nwithin AL",     age_CR   = "Age\nwithin CR",
                     interaction = "Interaction\n\u0394z(CR) - \u0394z(AL)")
gdot_size_range <- c(1.5, 5)   # point size, mm
gdot_row_mm     <- 11          # row height of the saved page; fixes the label offset
gdot_nudge      <- (max(gdot_size_range) / 2 + 1.8) / gdot_row_mm

gdot_sel <- unique(unlist(lapply(gdot_order, function(k) {
  r <- blk_rank[blk_rank$contrast == k, ]
  blk_pick_top(r$GO[order(-abs(r$change))], gdot_top_n)
})))

gdot_df <- blk_rank[blk_rank$contrast %in% gdot_order & blk_rank$GO %in% gdot_sel, ]
gdot_df$n_genes  <- unname(lengths(blk_sets)[gdot_df$GO])
gdot_df$contrast <- factor(gdot_df$contrast, levels = gdot_order,
                           labels = gdot_col_labels[gdot_order])

gdot_ord  <- names(sort(tapply(abs(gdot_df$change), gdot_df$GO, max)))  # ascending = bottom up
gdot_term <- go_term2name$TERM[match(gdot_ord, go_term2name$GO)]
gdot_term <- make.unique(gdot_term)
gdot_df$GO <- factor(gdot_df$GO, levels = gdot_ord, labels = gdot_term)

gdot_z_lim <- max(abs(gdot_df$change))
gdot_p_lim <- max(-log10(pmax(gdot_df$p_adj, 1e-10)))

p_blk_dot <- ggplot(gdot_df, aes(contrast, GO)) +
  geom_point(aes(size = -log10(pmax(p_adj, 1e-10)), fill = change),
             shape = 21, colour = "grey30", stroke = 0.3) +
  geom_text(aes(label = n_genes), nudge_y = gdot_nudge, size = 1.9,
            colour = "grey25") +
  scale_fill_gradient2(low = "#2166AC", mid = "grey93", high = "#B2182B",
                       midpoint = 0, limits = c(-gdot_z_lim, gdot_z_lim),
                       name = "\u0394z") +
  scale_size_continuous(range = gdot_size_range, limits = c(0, gdot_p_lim),
                        name = expression(-log[10]~"adj. p")) +
  scale_y_discrete(expand = expansion(add = c(0.6, 0.9)),
                   labels = scales::label_wrap(40)) +
  labs(x = NULL, y = NULL,
       title = sprintf("Top %d GO:BP terms per contrast by |\u0394z| in block connectivity",
                       gdot_top_n),
       subtitle = sprintf("Union of the per-contrast selections, every selected term shown in all %d contrasts.\nTerms ordered by the largest |\u0394z| across contrasts, numbers above the dots = genes in the term",
                          length(gdot_order)),
       caption = fig_caption(
         blk_celltype, method_labels[[reference_method]],
         "\u0394z = condition-mean z of the second group minus the first",
         sprintf("p = fraction of the %s size-matched decoy sets reaching |\u0394z| of the actual term in the\nsame contrast, BH adjusted within contrast (floor 1/%d)",
                 scales::comma(blk_n_bg), blk_n_bg + 1),
         sprintf("redundancy filter on each ranked list: overlap coefficient > %.2f or Jaccard > %.2f drops a term",
                 blk_max_overlap, blk_max_jaccard))) +
  theme(axis.text.y  = element_text(size = 6.5),
        axis.text.x  = element_text(size = 7, lineheight = 1.1),
        plot.caption = element_text(size = 6.5))

save_pdf(p_blk_dot, "G2_block_connectivity_dotplot.pdf",
         width = 8.6, height = nlevels(gdot_df$GO) * gdot_row_mm / 25.4 + 3.2,
         encoding = "Greek.enc")

# young-vs-old replication plot
blk_rep <- merge(blk_rank[blk_rank$contrast == "diet_young", c("GO", "change", "p")],
                 blk_rank[blk_rank$contrast == "diet_old",   c("GO", "change", "p")],
                 by = "GO", suffixes = c("_young", "_old"))
blk_rep <- blk_rep[blk_rep$p_young < 0.05 & blk_rep$p_old < 0.05 &
                     sign(blk_rep$change_young) == sign(blk_rep$change_old), ]

# strongest in both ages first, redundancy-filtered as in G1/G2
blk_rep_ord <- blk_rep$GO[order(-pmin(abs(blk_rep$change_young), abs(blk_rep$change_old)))]
blk_rep     <- blk_rep[match(blk_pick_top(blk_rep_ord, blk_rep_n_label), blk_rep$GO), ]
blk_rep$term <- sprintf("%s (%d)", go_term2name$TERM[match(blk_rep$GO, go_term2name$GO)],
                        lengths(blk_sets)[blk_rep$GO])
blk_rep$term <- factor(blk_rep$term, levels = rev(blk_rep$term))

blk_rep_long <- rbind(
  data.frame(term = blk_rep$term, age = "young", dz = blk_rep$change_young),
  data.frame(term = blk_rep$term, age = "old",   dz = blk_rep$change_old))
blk_rep_long$age <- factor(blk_rep_long$age, levels = c("young", "old"))

p_blk_rep <- ggplot(blk_rep, aes(y = term)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey55") +
  geom_segment(aes(x = change_young, xend = change_old, yend = term),
               colour = "grey60", linewidth = 0.6) +
  geom_point(data = blk_rep_long, aes(x = dz, fill = age), shape = 21, size = 2.8,
             colour = "grey30", stroke = 0.3) +
  scale_fill_manual(values = c(young = unname(cond_colours[["CR_young"]]),
                               old   = unname(cond_colours[["CR_old"]])),
                    name = "Diet contrast within") +
  scale_y_discrete(labels = scales::label_wrap(40)) +
  expand_limits(x = 0) +
  labs(x = "\u0394z (CR - AL)", y = NULL,
       title = "Diet effects reproduced in young and old animals",
       subtitle = sprintf("Top %d GO:BP terms at p < 0.05 in both contrasts, same sign, ordered by the smaller |\u0394z|, gene count in brackets",
                          blk_rep_n_label),
       caption = fig_caption(
         blk_celltype, method_labels[[reference_method]],
         sprintf("p = fraction of %s size-matched random sets reaching |\u0394z|",
                 scales::comma(blk_n_bg)),
         "redundancy filter as in G1/G2")) +
  theme(legend.position = "bottom",
        panel.grid.major.y = element_line(colour = "grey94", linewidth = 0.3))

save_pdf(p_blk_rep, "G3_diet_replication.pdf", width = 6.5, height = 4.5,
         encoding = "Greek.enc")

# ---- H PROFILE ENTROPY PER SAMPLE ----
ent_celltype <- blk_celltype
ent_gene_set <- NULL        # NULL = top3000 HVGs, conn_gene_set = top2000 HVGs
ent_top_k    <- 20          # strongest partners per gene entering the profile

ent_cms <- lapply(cor_mats[[reference_method]][[ent_celltype]][sample_levels],
                  function(cm) {
                    g <- if (is.null(ent_gene_set)) rownames(cm)
                    else intersect(ent_gene_set, rownames(cm))
                    cm[g, g, drop = FALSE]
                  })

# one partner set per gene from the mean |r| across all 8 samples: the selection
# is orthogonal to the contrasts and cannot carry a sample-specific noise floor
ent_keep <- local({
  a <- Reduce(`+`, lapply(ent_cms, abs)) / length(ent_cms)
  diag(a) <- 0
  a >= matrixStats::rowOrderStats(a, which = ncol(a) - ent_top_k + 1L)
})

ent_gene <- entropy_matrix(ent_cms, keep = ent_keep)

ent_df <- data.frame(sample = factor(sample_levels, levels = sample_levels),
                     H = colMeans(ent_gene)[sample_levels], row.names = NULL)
ent_df$condition <- factor(sample_to_condition[as.character(ent_df$sample)],
                           levels = condition_levels)
ent_df$age  <- factor(sub(".*_", "", as.character(ent_df$condition)),
                      levels = c("young", "old"))
ent_df$diet <- factor(sub("_.*", "", as.character(ent_df$condition)),
                      levels = c("AL", "CR"))

ent_cond <- dplyr::summarise(dplyr::group_by(ent_df, condition, age, diet),
                             H = mean(H), .groups = "drop")

p_ent <- ggplot(ent_df, aes(age, H)) +
  geom_line(data = ent_cond, aes(group = diet, colour = diet), linewidth = 0.6,
            show.legend = FALSE) +
  geom_point(aes(fill = condition), shape = 21, colour = "grey30", size = 2.6,
             stroke = 0.3,
             position = position_jitter(width = 0.07, height = 0, seed = seed)) +
  scale_colour_manual(values = blk_line_colours) +
  scale_cond_fill() +
  labs(x = NULL, y = sprintf("Shannon entropy (top %d partners)", ent_top_k),
       title = "Shannon entropy of the per-gene correlation matrices across conditions",
       subtitle = "Points = the 8 samples, coloured lines = condition means",
       caption = fig_caption(
         ent_celltype, method_labels[[reference_method]],
         sprintf("H = -sum(p*log(p)) over the %d strongest partners of a gene, p = |r| divided by their sum, diagonal excluded (%s genes)",
                 ent_top_k, scales::comma(nrow(ent_gene))),
         "\npartner set fixed per gene from the mean |r| across all 8 samples",
         "per sample: mean H over all genes")) +
  theme(legend.position = "bottom")

save_pdf(p_ent, "H1_profile_entropy.pdf", width = 8, height = 5.5)

# ---- I SEGMENT AXIS COMPRESSION ----
axis_gene_set  <- conn_gene_set          # same genes the F connectivity ran on
axis_celltype  <- paste0(pool_prefix_sc, "-raw")
axis_n_plot    <- 15000                  # cells per diet, plotting only
axis_n_label   <- 10                     # labelled genes in the loading panel

segment_colours <- setNames(c("#2166AC", "#E08214", "#1B7837"), base_celltypes)

# pooled z-scaling: one scale for all cells, so a per-sample variance difference
# is not normalised away before it can be measured
axis_genes <- intersect(axis_gene_set, rownames(seurat_obj))
X <- t(as.matrix(SeuratObject::LayerData(seurat_obj, assay = "RNA",
                                         layer = "data")[axis_genes, ]))
X <- scale(X, center = TRUE, scale = TRUE)
X <- X[, apply(X, 2, function(v) all(is.finite(v))), drop = FALSE]

axis_seg <- factor(as.character(seurat_obj$Annotation_lvl1), levels = base_celltypes)
axis_smp <- factor(as.character(seurat_obj$sample), levels = sample_levels)

# plane through the three segment centroids; diet enters nowhere
axis_ctr <- vapply(base_celltypes, function(s)
  colMeans(X[axis_seg == s, , drop = FALSE]), numeric(ncol(X)))
a1 <- axis_ctr[, "PT-S3"] - axis_ctr[, "PT-S1"]
a1 <- a1 / sqrt(sum(a1^2))
a2 <- axis_ctr[, "PT-S2"] - rowMeans(axis_ctr[, c("PT-S1", "PT-S3")])
a2 <- a2 - sum(a2 * a1) * a1                      # orthogonalise against a1
a2 <- a2 / sqrt(sum(a2^2))

add_diet_meta <- function(d) {
  d$sample    <- factor(as.character(d$sample), levels = sample_levels)
  d$condition <- factor(sample_to_condition[as.character(d$sample)],
                        levels = condition_levels)
  d$diet      <- factor(sub("_.*", "", as.character(d$condition)),
                        levels = c("AL", "CR"))
  d$age       <- factor(sub(".*_", "", as.character(d$condition)),
                        levels = c("young", "old"))
  d
}

axis_xy <- add_diet_meta(data.frame(
  x = as.numeric(X %*% a1), y = as.numeric(X %*% a2),
  segment = axis_seg, sample = axis_smp))
axis_xy$x <- axis_xy$x - mean(axis_xy$x)
axis_xy$y <- axis_xy$y - mean(axis_xy$y)

# distance between centroids, spread along the axis, spread in every other
# direction -- the three add up to the total within-sample variance
axis_decomp <- add_diet_meta(do.call(rbind, lapply(sample_levels, function(s) {
  i   <- which(axis_smp == s)
  g   <- droplevels(axis_seg[i])
  Y   <- X[i, , drop = FALSE]
  ctr <- rowsum(Y, g) / as.vector(table(g))
  res <- Y - ctr[as.character(g), , drop = FALSE]
  r1  <- as.numeric(res %*% a1)
  pc  <- setNames(as.numeric(ctr %*% a1), rownames(ctr))
  data.frame(sample = s,
             dist_axis   = pc[["PT-S3"]] - pc[["PT-S1"]],
             within_axis = mean(r1^2),
             within_orth = mean(rowSums(res^2)) - mean(r1^2),
             between     = mean(rowSums((ctr[as.character(g), , drop = FALSE] -
                                           rep(colMeans(Y), each = length(i)))^2)),
             row.names = NULL)
})))
axis_decomp$conn <- colMeans(conn_raw[[axis_celltype]])[as.character(axis_decomp$sample)]

axis_effect <- function(v) {
  m <- tapply(v, axis_decomp$diet, mean)
  100 * (m[["CR"]] / m[["AL"]] - 1)
}

# exact permutation over all balanced 4 vs 4 splits; smallest attainable p = 2/70
axis_perm_p <- function(v, diet) {
  obs    <- mean(v[diet == "CR"]) - mean(v[diet == "AL"])
  splits <- utils::combn(length(v), sum(diet == "CR"))
  null   <- apply(splits, 2, function(i) mean(v[i]) - mean(v[-i]))
  mean(abs(null) >= abs(obs) - 1e-12)
}

# dominant shared axis per sample, from the same correlation matrices as F
axis_g <- intersect(axis_gene_set, rownames(
  cor_mats[[reference_method]][[axis_celltype]][[sample_levels[1]]]))
axis_eig <- lapply(cor_mats[[reference_method]][[axis_celltype]][sample_levels],
                   function(cm) eigen(cm[axis_g, axis_g, drop = FALSE], symmetric = TRUE))
axis_L <- vapply(axis_eig, function(e) e$vectors[, 1], numeric(length(axis_g)))
rownames(axis_L) <- axis_g
colnames(axis_L) <- sample_levels
# eigenvector sign is arbitrary; anchor every sample on the first one
axis_L <- sweep(axis_L, 2, sign(stats::cor(axis_L, axis_L[, 1])[, 1]), "*")
axis_congruence <- stats::cor(axis_L)
stopifnot(all(axis_congruence > 0))      # fails if a sample carries another axis
axis_load <- rowMeans(axis_L)

axis_pc1 <- vapply(axis_eig, function(e) e$values[1] / sum(e$values), numeric(1))
axis_decomp$pc1 <- axis_pc1[as.character(axis_decomp$sample)]

# segment contrast of the same genes, computed without any correlation
axis_e   <- SeuratObject::LayerData(seurat_obj, assay = "RNA",
                                    layer = "data")[names(axis_load), ]
axis_lfc <- vapply(base_celltypes, function(s) {
  i <- axis_seg == s
  Matrix::rowMeans(axis_e[, i]) - Matrix::rowMeans(axis_e[, !i])
}, numeric(length(axis_load)))

axis_gene_df <- data.frame(gene = names(axis_load), load = as.numeric(axis_load),
                           s1s3 = as.numeric(axis_lfc[, "PT-S3"] -
                                               axis_lfc[, "PT-S1"]),
                           row.names = NULL)
axis_gene_r   <- stats::cor(axis_gene_df$load, axis_gene_df$s1s3)
axis_gene_rho <- stats::cor(axis_gene_df$load, axis_gene_df$s1s3, method = "spearman")

# --- I1 A: the shared axis is the segment gradient ---
p_axis_load <- ggplot(axis_gene_df, aes(s1s3, load)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey55") +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey55") +
  geom_point(size = 0.6, alpha = 0.4, colour = "grey35") +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, linewidth = 0.5,
              colour = "#B2182B") +
  ggrepel::geom_text_repel(
    data = axis_gene_df[order(-abs(axis_gene_df$load))[seq_len(axis_n_label)], ],
    aes(label = gene), size = 2.1, max.overlaps = Inf, seed = seed,
    segment.size = 0.2, colour = "grey20") +
  annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.6, size = 3, parse = TRUE,
           label = sprintf("italic(r) == %.2f*';'~~rho == %.2f", axis_gene_r, axis_gene_rho)) +
  labs(x = "Expression difference PT-S3 - PT-S1 (log-normalised)",
       y = "Loading on co-expression PC1",
       subtitle = "One point per gene. x from segment means only, y from the correlation matrices")

# --- I1 B: the plane, AL vs CR, AL outline repeated in the CR panel ---
axis_cen <- dplyr::summarise(dplyr::group_by(axis_xy, diet, segment),
                             x = mean(x), y = mean(y), .groups = "drop")
axis_cen <- axis_cen[order(axis_cen$diet, match(axis_cen$segment, base_celltypes)), ]

axis_plot_xy <- with_seed(seed, axis_xy[unlist(lapply(
  split(seq_len(nrow(axis_xy)), axis_xy$diet),
  function(i) sample(i, min(axis_n_plot, length(i))))), ])

axis_ref_xy  <- axis_plot_xy[axis_plot_xy$diet == "AL", ]
axis_ref_cen <- axis_cen[axis_cen$diet == "AL", ]
axis_ref_xy$diet  <- factor("CR", levels = levels(axis_plot_xy$diet))
axis_ref_cen$diet <- factor("CR", levels = levels(axis_cen$diet))

p_axis_plane <- ggplot(axis_plot_xy, aes(x, y)) +
  geom_point(aes(colour = segment), size = 0.25, alpha = 0.22) +
  stat_ellipse(data = axis_ref_xy, aes(group = segment, linetype = "AL outline"),
               level = 0.68, linewidth = 0.45, colour = "grey30") +
  geom_polygon(data = axis_ref_cen, fill = NA, colour = "grey30", linewidth = 0.45,
               linetype = "22", show.legend = FALSE) +
  stat_ellipse(aes(colour = segment), level = 0.68, linewidth = 0.5) +
  geom_polygon(data = axis_cen, aes(group = diet), fill = NA, colour = "grey15",
               linewidth = 0.45) +
  geom_point(data = axis_cen, aes(fill = segment), shape = 21, size = 2.6,
             colour = "grey15", stroke = 0.4) +
  facet_wrap(~ diet) +
  coord_equal() +
  scale_colour_manual(values = segment_colours, name = NULL) +
  scale_fill_manual(values = segment_colours, guide = "none") +
  scale_linetype_manual(values = c(`AL outline` = "22"), name = NULL) +
  guides(colour = guide_legend(override.aes = list(size = 2, alpha = 1))) +
  labs(x = "PT-S1 to PT-S3 axis (z units)", y = "PT-S2 axis (z units)",
       subtitle = sprintf("Both panels on one absolute scale. Ellipses = central 68%% of each segment, cells pooled over samples. %s cells per diet shown",
                          scales::comma(nrow(axis_plot_xy) / 2)))

# --- I1 C: the decomposition, per sample, relative to the AL mean ---
axis_comp_levels <- c(
  dist_axis   = "Distance between\nsegment centroids",
  within_axis = "Spread along\nthe axis",
  within_orth = "Spread in all\nother directions")

axis_rel <- add_diet_meta(do.call(rbind, lapply(names(axis_comp_levels), function(k) {
  v <- axis_decomp[[k]]
  data.frame(sample = axis_decomp$sample, component = k,
             rel = 100 * (v / mean(v[axis_decomp$diet == "AL"]) - 1), row.names = NULL)
})))
axis_rel$x <- match(axis_rel$component, names(axis_comp_levels)) +
  ifelse(axis_rel$diet == "AL", -0.2, 0.2)

axis_eff_lab <- data.frame(
  x   = seq_along(axis_comp_levels),
  lab = sprintf("%+.1f%%\np = %.3f",
                vapply(names(axis_comp_levels),
                       function(k) axis_effect(axis_decomp[[k]]), numeric(1)),
                vapply(names(axis_comp_levels),
                       function(k) axis_perm_p(axis_decomp[[k]], axis_decomp$diet), numeric(1))))
axis_diet_lab <- data.frame(
  x   = rep(seq_along(axis_comp_levels), each = 2) + c(-0.2, 0.2),
  lab = rep(c("AL", "CR"), length(axis_comp_levels)))

p_axis_comp <- ggplot(axis_rel, aes(x, rel)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey55") +
  stat_summary(aes(group = x), fun = mean, geom = "crossbar", width = 0.3,
               linewidth = 0.3, colour = "grey50") +
  geom_point(aes(fill = condition), shape = 21, colour = "grey30", size = 2.4,
             stroke = 0.3,
             position = position_jitter(width = 0.05, height = 0, seed = seed)) +
  geom_text(data = axis_eff_lab, aes(x, Inf, label = lab), inherit.aes = FALSE,
            vjust = 1.5, size = 3) +
  geom_text(data = axis_diet_lab, aes(x, -Inf, label = lab), inherit.aes = FALSE,
            vjust = -0.6, size = 2.6, colour = "grey35") +
  scale_x_continuous(breaks = seq_along(axis_comp_levels),
                     labels = unname(axis_comp_levels)) +
  scale_y_continuous(expand = expansion(mult = 0.15)) +
  scale_cond_fill() +
  labs(x = NULL, y = "Change relative to AL mean (%)",
       subtitle = "One point per sample, bar = diet mean, label = CR mean relative to AL mean") +
  theme(panel.grid.major.x = element_blank())

# --- I1 D: the compression is what the F connectivity was measuring ---
axis_link_r <- stats::cor(axis_decomp$dist_axis, axis_decomp$conn)

p_axis_link <- ggplot(axis_decomp, aes(dist_axis, conn)) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, linewidth = 0.5,
              colour = "grey55") +
  geom_point(aes(fill = condition), shape = 21, colour = "grey30", size = 3,
             stroke = 0.3) +
  ggrepel::geom_text_repel(aes(label = sample), size = 2.6, colour = "grey25",
                           seed = seed, show.legend = FALSE) +
  annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.6, size = 3,
           label = sprintf("r = %.2f, n = %d samples", axis_link_r, nrow(axis_decomp))) +
  scale_cond_fill() +
  labs(x = "PT-S1 to PT-S3 centroid distance (z units)",
       y = "Connectivity (mean |atanh(r)|)",
       subtitle = "One point per sample, connectivity as in section F") +
  guides(fill = "none")

save_pdf(patchwork::wrap_plots(
  patchwork::wrap_plots(p_axis_load, p_axis_plane, nrow = 1, widths = c(1, 1.6)),
  patchwork::wrap_plots(p_axis_comp, p_axis_link, nrow = 1, widths = c(1.4, 1)),
  ncol = 1, heights = c(1, 1), guides = "collect") +
    patchwork::plot_annotation(
      title = "PT-S1 to PT-S3 segment axis in AL and CR",
      subtitle = sprintf("%s, %s, %d genes of the top-%d HVG set. Axes defined from the segment centroids over all cells; diet does not enter their definition",
                         axis_celltype, method_labels[[reference_method]],
                         ncol(X), conn_n_hvg),
      caption = fig_caption(
        "genes z-scored once across all cells, not per sample",
        "axis 1 = S3 centroid - S1 centroid; axis 2 = S2 centroid - midpoint(S1, S3), orthogonalised against axis 1, both unit length",
        sprintf("loading = first eigenvector of the sample correlation matrix, sign-anchored on sample %s, averaged over samples (pairwise congruence %.2f-%.2f)",
                sample_levels[1], min(axis_congruence[lower.tri(axis_congruence)]),
                max(axis_congruence[lower.tri(axis_congruence)])),
        "spread along / orthogonal = mean squared within-segment residual projected on axis 1 and on its complement"),
      theme = theme_caption) &
    theme(legend.position = "bottom"),
  "I1_axis_compression.pdf", width = 12, height = 10)

# ---- RUN SUMMARY ----
pkg_version <- function(pkg)
  tryCatch(as.character(utils::packageVersion(pkg)), error = function(e) "not installed")

run_lines <- c(
  sprintf("data                     %s", data_path),
  sprintf("samples                  %s",
          paste(sprintf("%s (%s, %s)", sample_levels,
                        condition_labels[sample_to_condition[sample_levels]],
                        sample_to_batch[sample_levels]), collapse = ", ")),
  sprintf("cell types               %s", paste(base_celltypes, collapse = ", ")),
  "",
  sprintf("seed                     %d", seed),
  sprintf("gene selection           detection >= %.0f%%, top-%d HVGs per sample, voted within condition -> %d genes",
          100 * min_detection_fraction, n_hvg_per_sample, length(pt_gene_set)),
  sprintf("metacells                k = %d, %s | per sample %s",
          k_metacell, paste(metacell_methods, collapse = " / "),
          paste(sprintf("%s=%d", all_celltypes, n_eq_per_celltype), collapse = ", ")),
  sprintf("GO:BP term size          %d-%d (org.Mm.eg.db scale, ancestor-propagated)",
          go_min_term_size, go_max_term_size),
  sprintf("permutations per test    %s", scales::comma(n_perm)),
  "",
  sprintf("connectivity arms        %s x %s", paste(conn_celltypes, collapse = ", "),
          paste(names(conn_scalings), collapse = " / ")),
  "connectivity model       ~ age_c * diet_c, +-0.5 effect coding, eBayes(trend, robust), no batch term",
  sprintf("connectivity hits        %s", paste(sprintf("%s %s %d/%.0f",
                                                       conn_counts$arm, conn_counts$term, conn_counts$n_sig, conn_counts$exp_sig),
                                               collapse = " | ")),
  sprintf("set-level terms          %s of %s GO:BP terms with >= %d universe genes",
          scales::comma(length(blk_sets)), scales::comma(length(unique(go_term2gene$GO))),
          blk_min_genes),
  "",
  sprintf("output                   %s", normalizePath(out_dir, mustWork = FALSE)),
  sprintf("R %s | Seurat %s | hdWGCNA %s | finished %s",
          getRversion(), pkg_version("Seurat"), pkg_version("hdWGCNA"),
          format(Sys.time(), "%Y-%m-%d %H:%M")))

p_summary <- ggplot() +
  annotate("text", x = 0, y = 1, hjust = 0, vjust = 1, size = 5, fontface = "bold",
           label = "00_pipeline.R - run summary") +
  annotate("text", x = 0, y = 0.94, hjust = 0, vjust = 1, size = 3.1, family = "mono",
           label = paste(run_lines, collapse = "\n")) +
  scale_x_continuous(limits = c(0, 1)) +
  scale_y_continuous(limits = c(0, 1)) +
  theme_void()

save_pdf(p_summary, "00_run_summary.pdf", width = 18, height = 7)

# ---- THESIS MAIN FIGURES ----
thesis_main_dir <- file.path(out_dir, "thesis_figures")
dir.create(thesis_main_dir, showWarnings = FALSE)
stopifnot("CP1253.enc" %in% list.files(system.file("enc", package = "grDevices")))

thesis_w   <- 6.3      # A4 text width, inches
fig_sample <- "01"     # example sample in Fig. 6
sc_raw     <- paste0(pool_prefix_sc, "-raw")
focus_sets <- c(base_celltypes, sc_raw)
x_mid      <- (length(sample_levels) + 1) / 2   # centre of the sample axis

# printed sizes; CP1253 has no >=, x-sign or typographic minus
pt <- function(x) x / ggplot2::.pt
old_theme <- theme_set(theme_qc(base_size = 8) + theme(
  legend.position    = "bottom",
  legend.title       = element_text(size = 7),
  legend.text        = element_text(size = 7),
  legend.key.size    = unit(3, "mm"),
  legend.margin      = margin(0, 0, 0, 0),
  legend.box.spacing = unit(2, "mm"),
  strip.text         = element_text(size = 7, face = "bold", lineheight = 0.9),
  axis.text          = element_text(size = 6.5),
  axis.title         = element_text(size = 7.5),
  plot.tag           = element_text(size = 11, face = "bold")))

# two-line names where space is tight (x axes, 7 facets), one line otherwise
ds_lab  <- setNames(c(base_celltypes, "pool-\nmetacells-raw", "pool-\nmetacells-z",
                      "pool-\nsinglecells-raw", "pool-\nsinglecells-z"), all_celltypes)
ds_lab1 <- setNames(sub("\n", "", ds_lab, fixed = TRUE), all_celltypes)
mth_lab <- c(knn = "kNN metacells", random = "Random metacells")

relabel <- function(d, keep = all_celltypes, lab = ds_lab) {
  d <- d[as.character(d$celltype) %in% keep, , drop = FALSE]
  d$celltype <- factor(lab[as.character(d$celltype)], levels = unname(lab[keep]))
  if ("method" %in% names(d))
    d$method <- factor(mth_lab[as.character(d$method)], levels = mth_lab)
  d
}

# fewest decimals that keep the labels of one facet distinct
uniq_dec <- function(x) {
  d <- 2
  while (d < 4 && anyDuplicated(formatC(x[!is.na(x)], format = "f", digits = d))) d <- d + 1
  ifelse(is.na(x), NA_character_, formatC(x, format = "f", digits = d))
}

# first n GO terms in the given order whose gene sets are not redundant with a kept one
pick_nonredundant <- function(ids, sets, n) {
  keep <- character(0)
  for (id in ids) {
    if (length(keep) >= n) break
    if (!any(vapply(keep, function(k) blk_redundant(sets[[id]], sets[[k]]), logical(1))))
      keep <- c(keep, id)
  }
  keep
}

# horizontal colour bar under a panel
bar_below <- guide_colourbar(barwidth = unit(30, "mm"), barheight = unit(2.5, "mm"),
                             title.vjust = 1)

# titles and captions go into the Word legend; panel letters A, B, ...
save_thesis_fig <- function(p, file, height, width = thesis_w) {
  pages <- if (inherits(p, "gg")) list(p) else p
  clean <- labs(title = NULL, subtitle = NULL, caption = NULL)
  grDevices::pdf(file.path(thesis_main_dir, file), width = width, height = height,
                 onefile = TRUE, encoding = "CP1253.enc", useDingbats = FALSE)
  on.exit(grDevices::dev.off(), add = TRUE)
  for (pg in pages) {
    pg <- if (inherits(pg, "patchwork"))
      (pg & clean) + patchwork::plot_annotation(title = NULL, subtitle = NULL,
                                                caption = NULL, tag_levels = "A")
    else pg + clean
    print(pg)
  }
  invisible(NULL)
}

# --- shared panel builders ---
plot_go_stats <- function(dat, cnt) {
  relev <- function(d) {
    d$metric <- factor(as.integer(d$metric), levels = 1:2,
                       labels = c("Odds ratio\n(top 1% of |r|)", "AUC\n(all pairs)"))
    d
  }
  ggplot(relev(dat), aes(celltype, value, colour = condition)) +
    geom_hline(data = relev(go_stat_null), aes(yintercept = null),
               linetype = "dashed", colour = "grey55") +
    geom_point(position = position_jitter(width = 0.15, height = 0, seed = seed),
               size = 1.3, alpha = 0.9) +
    geom_text(data = relev(cnt), aes(celltype, Inf, label = label), inherit.aes = FALSE,
              vjust = 1.4, size = pt(6), colour = "grey30") +
    scale_y_continuous(labels = uniq_dec, expand = expansion(mult = c(0.05, 0.22))) +
    scale_cond_colour() +
    facet_grid(metric ~ method, scales = "free_y", switch = "y") +
    labs(x = NULL, y = NULL) +
    theme(strip.placement = "outside", axis.text.x = element_text(size = 6))
}

plot_rep <- function(d, label, breaks = rep_top_v) {
  p <- ggplot(d, aes(top_v, r, colour = condition, group = pair)) +
    geom_hline(yintercept = 0, linetype = "dotted", colour = "grey55") +
    geom_line(linewidth = 0.4, na.rm = TRUE) +
    geom_point(size = 0.9, na.rm = TRUE) +
    scale_x_log10(breaks = breaks,
                  expand = expansion(mult = c(0.05, if (label) 0.25 else 0.05))) +
    scale_cond_colour() +
    facet_grid(method ~ celltype) +
    labs(x = "Most variable genes considered (v, log scale)",
         y = "Correlation between replicates (r)") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  if (label)
    p <- p + ggrepel::geom_text_repel(
      data = d[d$top_v == max(rep_top_v) & !is.na(d$r), ], aes(label = n_genes),
      size = pt(5.5), direction = "y", nudge_x = 0.12, hjust = 0, segment.size = 0.15,
      min.segment.length = 0, box.padding = 0.1, seed = seed, show.legend = FALSE)
  p
}

term_short <- setNames(c("Age", "Diet", "Age x diet"), conn_term_labels)
plot_volcano <- function(ct, sc) {
  sel <- function(d) {
    d <- d[as.character(d$celltype) == ct & d$scaling == conn_scale_labels[[sc]], ]
    d$term <- factor(term_short[as.character(d$term)], levels = term_short)
    d
  }
  d <- sel(conn_stats)
  k <- sel(conn_counts)
  strip <- setNames(sprintf("%s\n%d / %d genes at p < %g",
                            k$term, k$n_sig, k$n_genes, conn_p),
                    as.character(k$term))
  top <- do.call(rbind, lapply(split(d, d$term), function(x) head(x[order(x$p), ], 4)))
  ggplot(d, aes(slope, -log10(p))) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey55") +
    geom_hline(yintercept = -log10(conn_p), linetype = "dotted", colour = "#B2182B",
               linewidth = 0.3) +
    geom_point(aes(colour = p < conn_p), size = 0.5, alpha = 0.6) +
    ggrepel::geom_text_repel(data = top, aes(label = gene), size = pt(5.5), seed = seed,
                             max.overlaps = Inf, segment.size = 0.15,
                             min.segment.length = 0, box.padding = 0.2) +
    facet_wrap(~ term, nrow = 1, scales = "free_y", labeller = as_labeller(strip)) +
    scale_y_continuous(expand = expansion(mult = c(0.03, 0.12))) +
    scale_colour_manual(values = c(`FALSE` = "grey70", `TRUE` = "#B2182B"), guide = "none") +
    labs(x = if (sc == "absolute") "Coefficient (change in connectivity)"
         else "Coefficient (change in connectivity relative to sample mean)",
         y = expression(-log[10]~p))
}

# sample-pair correlation heatmap of one dataset; text colour follows the printed value
plot_heat <- function(ct, row_labels = TRUE) {
  h   <- heatmap_df[as.character(heatmap_df$celltype) == ct, ]
  h$panel <- ds_lab1[[ct]]
  h$lv    <- round(100 * h$r) + 0
  mid <- mean(h$lv, na.rm = TRUE)
  ggplot(h, aes(x, y, fill = 100 * r)) +
    geom_tile(width = 1, height = 1, colour = "white", linewidth = 0.3) +
    geom_text(aes(label = ifelse(is.na(lv), "", sprintf("%.0f", lv)),
                  colour = !is.na(lv) & lv > mid),
              size = pt(6), show.legend = FALSE) +
    facet_wrap(~ panel) +
    scale_fill_gradient(low = "#F7FBFF", high = "#08306B", na.value = "grey92",
                        name = "r x 100", guide = bar_below) +
    scale_colour_manual(values = c(`FALSE` = "grey15", `TRUE` = "white"), guide = "none") +
    scale_x_continuous(breaks = axis_pos, labels = axis_label, expand = expansion(add = 0.6)) +
    scale_y_reverse(breaks = axis_pos, labels = if (row_labels) axis_label else NULL,
                    expand = expansion(add = 0.6)) +
    coord_fixed() +
    labs(x = NULL, y = NULL) +
    theme(panel.grid = element_blank(), axis.ticks = element_blank(),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, colour = axis_col,
                                     size = 6),
          axis.text.y = if (row_labels) element_text(colour = axis_col, size = 6)
          else element_blank())
}

# extra margin so ggrepel finds room for labels of points at the panel edge
pca_expand <- list(scale_x_continuous(expand = expansion(mult = 0.12)),
                   scale_y_continuous(expand = expansion(mult = 0.12)))
pca_panels <- function(ct) list(
  plot_sample_pca(sample_pca(metacells[[reference_method]][[ct]], "expression",
                             gene_sets[[ct]]), NULL) + pca_expand,
  plot_sample_pca(sample_pca(cor_mats[[reference_method]][[ct]], "structure"), NULL) +
    pca_expand)

# --- Fig. 3: gene selection ---
stage_short <- c("Expressed", "Detected (>= 10% of nuclei)",
                 "Variable (top 3,000 per sample)", "Final set (shared by segments)")
f3a_df <- gene_stage_df
f3a_df$stage <- factor(as.integer(f3a_df$stage), levels = 1:4, labels = stage_short)
f3a <- ggplot(f3a_df, aes(celltype, n, fill = stage)) +
  geom_col(position = position_dodge(0.85), width = 0.8) +
  geom_text(aes(label = scales::comma(n)), position = position_dodge(0.85),
            angle = 90, hjust = -0.1, size = pt(6)) +
  scale_fill_manual(values = setNames(c("grey80", "#9EC9C0", "#4E8FA6", "#1F3F5B"),
                                      stage_short), name = NULL) +
  scale_y_continuous(labels = scales::comma, expand = expansion(mult = c(0, 0.3))) +
  guides(fill = guide_legend(ncol = 1)) +
  labs(x = NULL, y = "Genes") +
  theme(legend.box.spacing = unit(0.5, "mm"))

f3b <- ggplot(hvg_curve_df[hvg_curve_df$var_std > 0, ],
              aes(rank, var_std, colour = condition, group = sample)) +
  geom_vline(xintercept = setdiff(hvg_sweep_n, n_hvg_per_sample), linetype = "dotted",
             colour = "grey65", linewidth = 0.3) +
  geom_vline(xintercept = n_hvg_per_sample, colour = "grey30", linewidth = 0.4) +
  geom_line(linewidth = 0.4) +
  facet_wrap(~ celltype, nrow = 1) +
  scale_x_continuous(breaks = c(0, 1000, 3000, 5000), labels = scales::comma) +
  scale_y_log10() +
  scale_cond_colour() +
  labs(x = "Gene rank within sample", y = "Standardised variance") +
  theme(legend.position = "none")

f3c_lab <- vapply(names(hvg_run_pairs), function(k) {
  p <- as.numeric(hvg_run_pairs[[k]])
  sprintf("%s vs %s HVGs\nslope %.2f", scales::comma(p[1]), scales::comma(p[2]),
          hvg_conn_fit$slope[as.character(hvg_conn_fit$pair) == k])
}, character(1))
f3c <- ggplot(hvg_conn_df, aes(conn_a, conn_b, colour = condition)) +
  geom_abline(linetype = "dashed", colour = "grey45") +
  geom_point(size = 0.25, alpha = 0.3) +
  facet_wrap(~ pair, nrow = 1, labeller = as_labeller(f3c_lab)) +
  scale_x_continuous(limits = hvg_conn_lim, breaks = seq(0.1, 0.4, 0.1)) +
  scale_y_continuous(limits = hvg_conn_lim, breaks = seq(0.1, 0.4, 0.1)) +
  coord_fixed() +
  scale_cond_colour() +
  guides(colour = guide_legend(override.aes = list(size = 2, alpha = 1))) +
  labs(x = "Connectivity, smaller cut-off", y = "Connectivity, larger cut-off")

# GO:BP over-representation of the selected genes under one universe: the n best
# non-redundant terms by FDR; order() is stable, so FDR ties keep enricher's raw p order
plot_go_selection <- function(universe, n = 10, fdr_lim = NULL, size_lim = NULL) {
  d    <- enrich_df[enrich_df$universe == universe, ]
  d    <- d[order(d$padj), ]
  t2g  <- go_term2gene[go_term2gene$GO %in% d$id, ]
  sets <- lapply(split(t2g$SYMBOL, t2g$GO), function(g) intersect(g, pt_gene_set))
  d    <- d[match(pick_nonredundant(d$id, sets, n), d$id), ]
  d$term <- factor(make.unique(d$term), levels = rev(make.unique(d$term)))
  ggplot(d, aes(fold, term)) +
    geom_vline(xintercept = 1, linetype = "dashed", colour = "grey55") +
    geom_point(aes(size = count, colour = -log10(padj))) +
    scale_y_discrete(labels = scales::label_wrap(30)) +
    scale_colour_gradient(low = "grey75", high = "#B2182B", limits = fdr_lim,
                          name = expression(-log[10]~FDR)) +
    scale_size_continuous(range = c(1, 3.5), limits = size_lim, name = "Genes") +
    labs(x = "Fold enrichment", y = NULL) +
    theme(legend.position = "right", axis.text.y = element_text(size = 6))
}

# main figure: universe of all genes expressed in any annotated kidney cell type
f3d <- plot_go_selection("detected")

save_thesis_fig(patchwork::wrap_plots(A = f3a, B = f3b, C = f3c, D = patchwork::free(f3d),
                                      design = "AD\nBB\nCC", widths = c(0.85, 1.15),
                                      heights = c(1.3, 0.75, 0.8)),
                "Fig03_gene_selection.pdf", height = 8.4)

# --- Fig. 5: cells, metacells, correlation distributions, segment purity ---
f5a_lab <- data.frame(
  celltype = factor(base_celltypes, levels = base_celltypes),
  lab = sprintf("%d kept per sample",
                tapply(cell_count_df$count_after, cell_count_df$celltype, min)[base_celltypes]))
f5a <- ggplot(cell_count_df, aes(sample)) +
  geom_col(aes(y = count_before), fill = "grey85", width = 0.75) +
  geom_col(aes(y = count_after, fill = condition), width = 0.75) +
  geom_text(data = f5a_lab, aes(x_mid, Inf, label = lab), inherit.aes = FALSE,
            hjust = 0.5, vjust = 1.5, size = pt(6.5)) +
  facet_wrap(~ celltype, nrow = 1) +
  scale_cond_fill() +
  scale_y_continuous(labels = scales::comma, expand = expansion(mult = c(0, 0.12))) +
  labs(x = NULL, y = "Cells") +
  theme(axis.text.x = element_text(size = 6))

f5b_df <- relabel(mc_count_df, focus_sets, ds_lab1)
f5b_df$sample <- factor(as.character(f5b_df$sample), levels = sample_levels)
f5b <- ggplot(f5b_df, aes(sample)) +
  geom_col(aes(y = n_possible), fill = "grey88", width = 0.8) +
  geom_col(aes(y = n_built), fill = "grey65", width = 0.8) +
  geom_col(aes(y = n_kept, fill = condition), width = 0.8) +
  geom_text(data = unique(f5b_df[, c("celltype", "method", "n_kept")]),
            aes(x_mid, Inf, label = sprintf("%d kept per sample", n_kept)),
            inherit.aes = FALSE, hjust = 0.5, vjust = 1.5, size = pt(6.5)) +
  facet_grid(method ~ celltype) +
  scale_cond_fill(guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(x = NULL, y = "Metacells per sample") +
  theme(axis.text.x = element_text(size = 5.5, angle = 90, vjust = 0.5))

f5c <- ggplot(relabel(cor_dist_df[cor_dist_df$method == "knn", ], focus_sets, ds_lab1),
              aes(abs_cor, colour = condition, group = sample)) +
  geom_density(linewidth = 0.4, adjust = 1.5) +
  facet_wrap(~ celltype, nrow = 1, scales = "free_y") +
  coord_cartesian(xlim = c(0, 0.75)) +
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75)) +
  scale_cond_colour(guide = "none") +
  labs(x = "|r| between gene pairs (kNN metacells)", y = "Density")

f5d_df <- comp_df[as.character(comp_df$celltype) == sc_raw, ]
f5d_df$method <- factor(mth_lab[as.character(f5d_df$method)], levels = mth_lab)
f5d_df$sample <- factor(as.character(f5d_df$sample), levels = sample_levels)
f5d <- ggplot(f5d_df, aes(sample, purity, fill = condition)) +
  geom_violin(scale = "width", linewidth = 0.2, colour = "grey45") +
  geom_hline(yintercept = purity_baseline, linetype = "dashed", colour = "grey40") +
  facet_wrap(~ method, ncol = 1) +
  scale_cond_fill(guide = "none") +
  scale_y_continuous(limits = c(NA, 1)) +
  labs(x = NULL, y = "Largest segment fraction per metacell") +
  theme(axis.text.x = element_text(size = 6))

save_thesis_fig(patchwork::wrap_plots(A = f5a, B = patchwork::free(f5b),
                                      C = patchwork::free(f5c), D = f5d,
                                      design = "AD\nBB\nCC", widths = c(1.7, 1),
                                      heights = c(1.1, 1.3, 0.75), guides = "collect") &
                  theme(legend.position = "bottom"),
                "Fig05_metacells.pdf", height = 8.4)

# --- Fig. 6: functional coherence ---
f6_edge <- go_edge_df[go_edge_df$method == reference_method &
                        as.character(go_edge_df$celltype) %in% focus_sets &
                        go_edge_df$sample == fig_sample, ]
f6a_df <- do.call(rbind, lapply(seq_len(nrow(f6_edge)), function(i) {
  s  <- f6_edge[i, ]
  n  <- as.double(s$a + s$b + s$c + s$d)
  nh <- as.double(s$n_high); ns <- as.double(s$n_shared)
  count    <- c(s$a, s$b, s$c, s$d)
  expected <- c(nh * ns, nh * (n - ns), (n - nh) * ns, (n - nh) * (n - ns)) / n
  data.frame(celltype = as.character(s$celltype),
             go  = factor(c("shared", "none", "shared", "none"), levels = c("shared", "none")),
             cor = factor(c("top 1%", "top 1%", "rest", "rest"), levels = c("rest", "top 1%")),
             ratio = log2(pmax(count, 0.5) / expected),
             lab = scales::comma(count),
             stringsAsFactors = FALSE)
}))
f6a_df <- relabel(f6a_df, focus_sets, ds_lab1)
f6a_strip <- setNames(
  sprintf("%s\nOR = %.2f, p = %s", ds_lab1[as.character(f6_edge$celltype)],
          f6_edge$odds_ratio,
          format.pval(f6_edge$perm_p, digits = 2, eps = 1 / (n_perm + 1))),
  ds_lab1[as.character(f6_edge$celltype)])
f6a_lim <- max(abs(f6a_df$ratio))
f6a <- ggplot(f6a_df, aes(go, cor, fill = ratio)) +
  geom_tile(colour = "white", linewidth = 0.6) +
  geom_text(aes(label = lab), size = pt(6)) +
  facet_wrap(~ celltype, nrow = 1, labeller = as_labeller(f6a_strip)) +
  scale_fill_gradient2(low = "#2166AC", mid = "grey95", high = "#B2182B",
                       limits = c(-f6a_lim, f6a_lim), name = "log2 observed / expected",
                       guide = bar_below) +
  scale_x_discrete(position = "top", expand = c(0, 0),
                   labels = c(shared = "shared\nGO term", none = "no shared\nterm")) +
  scale_y_discrete(expand = c(0, 0)) +
  labs(x = NULL, y = "|r| class") +
  theme(panel.grid = element_blank(), strip.placement = "outside",
        axis.ticks = element_blank())

f6b_df <- do.call(rbind, lapply(focus_sets, function(ct) {
  pv <- go_pair_vectors(cor_mats[[reference_method]][[ct]][[fig_sample]], go_annot)
  do.call(rbind, lapply(c(TRUE, FALSE), function(is_shared) {
    v <- pv$abs_cor[pv$shared == is_shared]
    data.frame(celltype = ct,
               class = if (is_shared) "shared GO:BP term" else "no shared term",
               abs_cor = with_seed(seed, sample(v, min(n_dist_subsample, length(v)))),
               stringsAsFactors = FALSE)
  }))
}))
f6b <- ggplot(relabel(f6b_df, focus_sets, ds_lab1), aes(abs_cor, colour = class)) +
  geom_density(linewidth = 0.5, adjust = 1.3) +
  facet_wrap(~ celltype, nrow = 1, scales = "free_y") +
  coord_cartesian(xlim = c(0, 0.75)) +
  scale_colour_manual(values = c("shared GO:BP term" = "#B2182B",
                                 "no shared term" = "#4393C3"), name = NULL) +
  labs(x = sprintf("|r| between gene pairs (sample %s)", fig_sample), y = "Density")

f6c <- ggplot(relabel(go_sweep_df[go_sweep_df$method == "knn", ], focus_sets, ds_lab1),
              aes(top_q, odds_ratio)) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey55") +
  geom_vline(xintercept = go_top_edge_q, linetype = "dotted", colour = "grey30") +
  geom_line(aes(colour = condition, group = sample), linewidth = 0.3, alpha = 0.6) +
  geom_line(data = relabel(sweep_median[sweep_median$method == "knn", ], focus_sets, ds_lab1),
            colour = "black", linewidth = 0.6) +
  facet_wrap(~ celltype, nrow = 1) +
  scale_x_log10(breaks = c(0.01, 0.1, 0.5), labels = c("1%", "10%", "50%")) +
  scale_cond_colour() +
  labs(x = "Fraction of gene pairs classified as highly correlated (log scale)",
       y = "Odds ratio") +
  theme(legend.position = "none", panel.spacing.x = unit(3, "mm"))

f6d <- plot_go_stats(relabel(go_stat_df, focus_sets), relabel(go_stat_counts, focus_sets))

save_thesis_fig(patchwork::wrap_plots(A = f6a, B = f6b, C = f6c, D = f6d,
                                      design = "A\nB\nC\nD",
                                      heights = c(0.95, 0.75, 0.75, 1.5)),
                "Fig06_GO_coherence.pdf", height = 8.6)

# --- Fig. 7: reproducibility ---
save_thesis_fig(patchwork::wrap_plots(
  A = plot_rep(relabel(rep_cor_df, focus_sets, ds_lab1), TRUE),
  B = plot_heat("PT-S1"), C = plot_heat(sc_raw, row_labels = FALSE),
  design = "AA\nBC", heights = c(1.25, 1)),
  "Fig07_reproducibility.pdf", height = 8.2)

# --- Fig. 8: sample-level PCA ---
save_thesis_fig(patchwork::wrap_plots(pca_panels(sc_raw), nrow = 1, guides = "collect") &
                  theme(legend.position = "bottom", legend.spacing.x = unit(8, "mm")),
                "Fig08_sample_PCA.pdf", height = 3.7)

# --- Fig. 9: connectivity model ---
# condition legend comes from panel B only
f9a <- ggplot(ent_df, aes(age, H)) +
  geom_line(data = ent_cond, aes(group = diet, colour = diet), linewidth = 0.5) +
  geom_point(aes(fill = condition), shape = 21, colour = "grey30", size = 2,
             stroke = 0.3, position = position_jitter(width = 0.08, height = 0, seed = seed)) +
  scale_colour_manual(values = blk_line_colours, guide = "none") +
  scale_cond_fill(guide = "none") +
  labs(x = NULL, y = "Mean profile entropy")

f9_keys  <- c("expr", "variance", "conn")
f9_names <- c("Expression", "Variance", "Connectivity")
f9b_df <- mol_dist[mol_dist$metric %in% mol_metric_levels[f9_keys], ]
f9b_df$metric <- factor(f9_names[match(as.character(f9b_df$metric),
                                       mol_metric_levels[f9_keys])], levels = f9_names)
f9b_mean <- aggregate(value ~ metric + sample, data = f9b_df, FUN = mean)
f9b <- ggplot(f9b_df, aes(sample, value)) +
  geom_violin(aes(fill = condition), scale = "width", linewidth = 0.2, colour = "grey45") +
  geom_errorbar(data = f9b_mean, aes(x = sample, ymin = value, ymax = value),
                inherit.aes = FALSE, width = 0.7, linewidth = 0.35) +
  facet_wrap(~ metric, nrow = 1, scales = "free_y") +
  scale_cond_fill() +
  labs(x = NULL, y = "Value per gene") +
  theme(axis.text.x = element_text(size = 6, angle = 90, vjust = 0.5))

gsea_sets_conn <- local({
  hit <- go_term2gene$SYMBOL %in% conn_gene_set
  lapply(split(go_term2gene$SYMBOL[hit], go_term2gene$GO[hit]), unique)
})
lol_df <- gsea_df[gsea_df$ct == sc_raw & gsea_df$rank == "t", ]
lol_df <- do.call(rbind, lapply(split(lol_df, list(lol_df$sc, lol_df$term), drop = TRUE),
                                function(d) {
                                  d <- d[order(d$p), ]
                                  d[match(pick_nonredundant(d$ID, gsea_sets_conn, 5), d$ID), ]
                                }))
lol_panels <- as.vector(t(outer(gsea_term_short[conn_term_levels], names(conn_scalings),
                                function(t, s) sprintf("%s | %s", t, s))))
lol_df$panel <- factor(sprintf("%s | %s", gsea_term_short[lol_df$term], lol_df$sc),
                       levels = lol_panels)
lol_df$key <- paste(lol_df$panel, lol_df$ID)
lol_df$key <- factor(lol_df$key, levels = lol_df$key[order(lol_df$NES)])
f9e <- ggplot(lol_df, aes(NES, key)) +
  geom_vline(xintercept = 0, colour = "grey55") +
  geom_segment(aes(x = 0, xend = NES, yend = key), colour = "grey60", linewidth = 0.5) +
  geom_point(aes(size = -log10(p), fill = NES), shape = 21, colour = "grey30", stroke = 0.3) +
  facet_wrap(~ panel, ncol = 2, scales = "free_y") +
  scale_y_discrete(labels = setNames(scales::label_wrap(60)(lol_df$Description),
                                     as.character(lol_df$key))) +
  scale_fill_gradient2(low = "#2166AC", mid = "grey93", high = "#B2182B",
                       midpoint = 0, guide = "none") +
  scale_size_continuous(range = c(1, 3), name = expression(-log[10]~p~"(nominal)")) +
  labs(x = "Normalised enrichment score (genes ranked by moderated t)", y = NULL) +
  theme(axis.text.y = element_text(size = 6, lineheight = 0.85),
        legend.location = "plot")

f9_top <- patchwork::wrap_plots(f9a, f9b, nrow = 1, widths = c(0.35, 1),
                                guides = "collect") &
  theme(legend.position = "bottom")
save_thesis_fig(patchwork::wrap_plots(
  f9_top, plot_volcano(sc_raw, "absolute"), plot_volcano(sc_raw, "relative"),
  patchwork::free(f9e), ncol = 1, heights = c(0.85, 0.8, 0.8, 2)),
  "Fig09_connectivity_model.pdf", height = 9)

# --- Fig. 10: term-level block connectivity ---
f10a <- ggplot(blk_rep, aes(y = term)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey55") +
  geom_segment(aes(x = change_young, xend = change_old, yend = term),
               colour = "grey60", linewidth = 0.5) +
  geom_point(data = blk_rep_long, aes(x = dz, fill = age), shape = 21, colour = "grey30",
             size = 2.2, stroke = 0.3) +
  scale_fill_manual(values = c(young = unname(cond_colours[["CR_young"]]),
                               old   = unname(cond_colours[["CR_old"]])),
                    name = "Diet contrast within") +
  scale_y_discrete(labels = scales::label_wrap(45)) +
  expand_limits(x = 0) +
  labs(x = "\u0394z (CR - AL)", y = NULL)

# examples: four of the ten replicated terms, chosen to show different but plausible
# patterns of change (ranks 1, 3, 7, 10 by min(dz young, dz old))
fig10_terms <- c("response to fatty acid",
                 "positive regulation of actin filament bundle assembly",
                 "monoatomic anion transmembrane transport",
                 "glucose metabolic process")
fig10_go <- go_term2name$GO[match(fig10_terms, go_term2name$TERM)]
stopifnot(!anyNA(fig10_go), all(fig10_go %in% names(blk_sets)))
as_fig10 <- function(d) {
  d <- d[d$GO %in% fig10_go, ]
  d$term <- factor(go_term2name$TERM[match(d$GO, go_term2name$GO)], levels = fig10_terms)
  d
}
f10b_cond <- as_fig10(blk_cond)
f10b <- ggplot(as_fig10(blk_df), aes(age, obs)) +
  geom_line(data = f10b_cond, aes(y = bg_mean, group = diet), colour = "grey70",
            linetype = "22", linewidth = 0.4) +
  geom_line(data = f10b_cond, aes(group = diet, colour = diet), linewidth = 0.5) +
  geom_point(aes(fill = condition), shape = 21, colour = "grey30", size = 1.8, stroke = 0.3,
             position = position_jitter(width = 0.07, height = 0, seed = seed)) +
  facet_wrap(~ term, nrow = 1, scales = "free_y", labeller = label_wrap_gen(20)) +
  scale_colour_manual(values = blk_line_colours, name = "Diet (condition mean)") +
  scale_cond_fill() +
  labs(x = NULL, y = "Block connectivity\n(mean |atanh(r)|)")

save_thesis_fig(patchwork::wrap_plots(A = f10a, B = patchwork::free(f10b),
                                      design = "A\nB", heights = c(1.15, 0.85)),
                "Fig10_term_connectivity.pdf", height = 7)

# --- Fig. 11: segmental axis ---
f11b <- p_axis_plane +
  facet_wrap(~ diet, ncol = 1) +
  guides(colour = guide_legend(nrow = 1, override.aes = list(size = 2, alpha = 1)),
         linetype = guide_legend(nrow = 1)) +
  theme(legend.box = "horizontal", legend.spacing.x = unit(1, "mm"),
        legend.key.spacing.x = unit(0.5, "mm"))
f11c <- p_axis_comp +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.45))) +
  guides(fill = guide_legend(nrow = 1))
f11_bottom <- patchwork::wrap_plots(f11c, p_axis_link, nrow = 1, widths = c(1.35, 1),
                                    guides = "collect") &
  theme(legend.position = "bottom")
save_thesis_fig(patchwork::wrap_plots(
  patchwork::wrap_plots(p_axis_load, f11b, nrow = 1, widths = c(1.3, 1)),
  f11_bottom, ncol = 1, heights = c(1.35, 1)),
  "Fig11_axis_compression.pdf", height = 8.4)

# --- supplementary figures ---
sc_pools <- grep(pool_prefix_sc, all_celltypes, value = TRUE)
s1c_df <- relabel(comp_df, sc_pools, ds_lab1)
s1c_df$sample <- factor(as.character(s1c_df$sample), levels = sample_levels)
s1c <- ggplot(s1c_df, aes(sample, purity, fill = condition)) +
  geom_violin(scale = "width", linewidth = 0.2, colour = "grey45") +
  geom_hline(yintercept = purity_baseline, linetype = "dashed", colour = "grey40") +
  facet_grid(method ~ celltype) +
  scale_cond_fill(guide = "none") +
  scale_y_continuous(limits = c(NA, 1)) +
  labs(x = NULL, y = "Largest segment\nfraction per metacell") +
  theme(axis.text.x = element_text(size = 6), strip.text = element_text(size = 6))
save_thesis_fig(patchwork::wrap_plots(
  A = plot_go_stats(relabel(go_stat_df), relabel(go_stat_counts)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 5.5)),
  B = plot_rep(relabel(rep_cor_df), FALSE, breaks = c(20, 100, 1000)) +
    guides(colour = "none") +
    theme(strip.text = element_text(size = 6)),
  C = s1c, design = "A\nB\nC", heights = c(1.35, 1.2, 0.85), guides = "collect") &
    theme(legend.position = "bottom"),
  "SFig01_pooling_comparison.pdf", height = 9)

save_thesis_fig(patchwork::wrap_plots(unlist(lapply(base_celltypes, pca_panels),
                                             recursive = FALSE),
                                      ncol = 2, guides = "collect") &
                  theme(legend.position = "bottom", legend.spacing.x = unit(8, "mm")),
                "SFig02_sample_PCA_segments.pdf", height = 8.6)

s3_keys  <- c("lib_pre_all", "lib_pre_pt", "lib_post_all", "lib_post_pt")
s3_names <- c("Before depth norm.,\nall cell types", "Before depth norm.,\nPT (equalised)",
              "After depth norm.,\nall cell types", "After depth norm.,\nPT (equalised)")
s3_df <- mol_dist[mol_dist$metric %in% mol_metric_levels[s3_keys], ]
s3_df$metric <- factor(s3_names[match(as.character(s3_df$metric), mol_metric_levels[s3_keys])],
                       levels = s3_names)
s3_x <- scale_x_discrete(labels = function(s) paste0(s, " (", sample_to_batch[s], ")"))
s3_theme <- theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 5.5))
s3a <- ggplot(s3_df, aes(sample, value)) +
  geom_violin(aes(fill = condition), scale = "width", linewidth = 0.2, colour = "grey45") +
  stat_summary(fun = median, geom = "crossbar", width = 0.6, linewidth = 0.25) +
  facet_wrap(~ metric, nrow = 1) +
  s3_x + scale_cond_fill() +
  labs(x = NULL, y = "Library size (log10 UMIs)") + s3_theme
s3b <- ggplot(mol_stats, aes(sample, singleton, fill = condition)) +
  geom_col(width = 0.7) +
  s3_x + scale_cond_fill(guide = "none") +
  labs(x = NULL, y = "Fraction of molecules\nwith one read") + s3_theme
save_thesis_fig(patchwork::wrap_plots(s3a, s3b, nrow = 1, widths = c(4, 1),
                                      guides = "collect") &
                  theme(legend.position = "bottom"),
                "SFig03_sequencing_depth.pdf", height = 3.6)

mc_raw <- paste0(pool_prefix_mc, "-raw")
save_thesis_fig(patchwork::wrap_plots(plot_volcano(mc_raw, "absolute"),
                                      plot_volcano(mc_raw, "relative"), ncol = 1),
                "SFig04_connectivity_model_metacell_pool.pdf", height = 5)

save_thesis_fig(lapply(names(blk_contrasts), blk_page),
                "SFig05_block_connectivity.pdf", width = 9.4, height = 6)

save_thesis_fig(p_blk_dot, "SFig06_block_connectivity_overview.pdf", height = 8.5)

theme_set(old_theme)