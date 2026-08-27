args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 7) {
  stop("usage: dada2_profile.R LIBRARY_ID ASSAY LAYOUT R1 R2 TAXONOMY OUTDIR")
}

library_id <- args[[1]]
assay <- args[[2]]
layout <- args[[3]]
r1 <- args[[4]]
r2 <- args[[5]]
taxonomy_reference <- args[[6]]
outdir <- args[[7]]

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
input_connection <- gzfile(r1, open = "rt")
first_record <- readLines(input_connection, n = 1, warn = FALSE)
close(input_connection)
if (length(first_record) == 0) {
  write.table(
    data.frame(
      library_id = character(), assay = character(), feature_id = character(),
      original_name = character(), original_rank = character(), abundance = numeric()
    ),
    file.path(outdir, "profile.tsv"), sep = "\t", quote = FALSE, row.names = FALSE
  )
  write.table(
    data.frame(feature_id = character(), count = numeric()),
    file.path(outdir, "asv_counts.tsv"), sep = "\t", quote = FALSE, row.names = FALSE
  )
  file.create(file.path(outdir, "asvs.fasta"))
  quit(save = "no", status = 0)
}

if (!requireNamespace("dada2", quietly = TRUE)) {
  stop("The R package 'dada2' is required")
}

valid_error_matrix <- function(matrix) {
  is.matrix(matrix) && nrow(matrix) == 16 && ncol(matrix) > 0 &&
    all(is.finite(matrix)) && all(matrix >= 0) && all(matrix <= 1)
}

learn_error_model <- function(path, direction) {
  messages <- character()
  learned <- withCallingHandlers(
    tryCatch(
      dada2::learnErrors(path, multithread = TRUE),
      error = function(error) {
        messages <<- c(messages, conditionMessage(error))
        NULL
      }
    ),
    warning = function(warning) {
      messages <<- c(messages, conditionMessage(warning))
      invokeRestart("muffleWarning")
    }
  )
  if (valid_error_matrix(learned)) {
    return(list(
      error_matrix = learned,
      method = "learned_per_library",
      reason = ""
    ))
  }

  fallback <- if (direction == "forward") {
    dada2::errBalancedF
  } else {
    dada2::errBalancedR
  }
  if (!valid_error_matrix(fallback)) {
    stop("DADA2 supplied neither a learned nor a valid bundled error matrix")
  }
  reason <- paste(unique(messages[nzchar(messages)]), collapse = "; ")
  if (!nzchar(reason)) {
    reason <- "learnErrors returned no valid error matrix"
  }
  warning(
    sprintf(
      "%s reads could not support a per-library error fit; using DADA2's bundled MiSeq 2x250 empirical %s matrix (%s)",
      direction, direction, reason
    ),
    call. = FALSE
  )
  list(
    error_matrix = fallback,
    method = paste0("dada2_errBalanced", if (direction == "forward") "F" else "R"),
    reason = reason
  )
}

filtered_r1 <- file.path(outdir, "filtered_R1.fastq.gz")
filtered_r2 <- file.path(outdir, "filtered_R2.fastq.gz")

if (layout == "paired") {
  dada2::filterAndTrim(
    r1, filtered_r1, r2, filtered_r2,
    maxN = 0, maxEE = c(2, 2), truncQ = 2,
    rm.phix = TRUE, compress = TRUE, multithread = TRUE
  )
  model_r1 <- learn_error_model(filtered_r1, "forward")
  model_r2 <- learn_error_model(filtered_r2, "reverse")
  derep_r1 <- dada2::derepFastq(filtered_r1)
  derep_r2 <- dada2::derepFastq(filtered_r2)
  denoised_r1 <- dada2::dada(
    derep_r1, err = model_r1$error_matrix, multithread = TRUE
  )
  denoised_r2 <- dada2::dada(
    derep_r2, err = model_r2$error_matrix, multithread = TRUE
  )
  merged <- dada2::mergePairs(denoised_r1, derep_r1, denoised_r2, derep_r2)
  table <- dada2::makeSequenceTable(merged)
  error_models <- data.frame(
    read = c("R1", "R2"),
    method = c(model_r1$method, model_r2$method),
    reason = c(model_r1$reason, model_r2$reason),
    stringsAsFactors = FALSE
  )
} else {
  dada2::filterAndTrim(
    r1, filtered_r1, maxN = 0, maxEE = 2, truncQ = 2,
    rm.phix = TRUE, compress = TRUE, multithread = TRUE
  )
  model_r1 <- learn_error_model(filtered_r1, "forward")
  derep <- dada2::derepFastq(filtered_r1)
  denoised <- dada2::dada(
    derep, err = model_r1$error_matrix, multithread = TRUE
  )
  table <- dada2::makeSequenceTable(denoised)
  error_models <- data.frame(
    read = "R1",
    method = model_r1$method,
    reason = model_r1$reason,
    stringsAsFactors = FALSE
  )
}

write.table(
  error_models, file.path(outdir, "error_model.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

table_nochim <- dada2::removeBimeraDenovo(
  table, method = "consensus", multithread = TRUE
)
taxa <- dada2::assignTaxonomy(
  table_nochim, taxonomy_reference, multithread = TRUE
)

sequences <- colnames(table_nochim)
feature_ids <- sprintf("ASV%06d", seq_along(sequences))
counts <- as.numeric(table_nochim[1, ])
rank_columns <- colnames(taxa)

deepest_name <- function(row) {
  values <- row[!is.na(row) & nzchar(row)]
  if (length(values) == 0) "Unclassified" else tail(values, 1)
}
deepest_rank <- function(row) {
  present <- which(!is.na(row) & nzchar(row))
  if (length(present) == 0) "unclassified" else tolower(rank_columns[tail(present, 1)])
}

profile <- data.frame(
  library_id = library_id,
  assay = assay,
  feature_id = feature_ids,
  original_name = apply(taxa, 1, deepest_name),
  original_rank = apply(taxa, 1, deepest_rank),
  abundance = counts,
  stringsAsFactors = FALSE
)
write.table(
  profile, file.path(outdir, "profile.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
write.table(
  data.frame(feature_id = feature_ids, count = counts),
  file.path(outdir, "asv_counts.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
fasta <- unlist(Map(
  function(id, sequence) c(paste0(">", id), sequence),
  feature_ids, sequences
))
writeLines(fasta, file.path(outdir, "asvs.fasta"))
