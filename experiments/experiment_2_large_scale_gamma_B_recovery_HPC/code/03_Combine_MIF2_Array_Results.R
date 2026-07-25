# ============================================================
# Combine results from the nine MIF2 array tasks
#
# Outputs:
#   results/combined_mif2_results.csv
#   results/best_fit.csv
#   results/best_mif2_object.rds
# ============================================================

options(stringsAsFactors = FALSE)

results_root <- file.path(
  "results",
  "array_output"
)

result_files <- list.files(
  path = results_root,
  pattern = "mif2_result\\.csv$",
  recursive = TRUE,
  full.names = TRUE
)

if (length(result_files) == 0) {
  stop(
    paste(
      "No task result files were found under",
      results_root
    )
  )
}

if (length(result_files) != 9) {
  warning(
    "Expected 9 task result files, but found ",
    length(result_files),
    "."
  )
}

mif_results <- do.call(
  rbind,
  lapply(
    result_files,
    read.csv
  )
)

rownames(mif_results) <- NULL

mif_results <- mif_results[
  order(
    mif_results$task_id
  ),
  ,
  drop = FALSE
]

dir.create(
  "results",
  recursive = TRUE,
  showWarnings = FALSE
)

combined_file <- file.path(
  "results",
  "combined_mif2_results.csv"
)

write.csv(
  mif_results,
  combined_file,
  row.names = FALSE
)

valid_results <- mif_results[
  is.finite(
    mif_results$logLik
  ),
  ,
  drop = FALSE
]

if (nrow(valid_results) == 0) {
  stop(
    "No task has a finite evaluated log-likelihood."
  )
}

best_fit <- valid_results[
  which.max(
    valid_results$logLik
  ),
  ,
  drop = FALSE
]

best_fit_file <- file.path(
  "results",
  "best_fit.csv"
)

write.csv(
  best_fit,
  best_fit_file,
  row.names = FALSE
)

best_task_id <- best_fit$task_id[[1]]

best_object_source <- file.path(
  results_root,
  sprintf(
    "task_%03d",
    best_task_id
  ),
  "mif2_object.rds"
)

if (!file.exists(best_object_source)) {
  stop(
    paste(
      "Best MIF2 object not found:",
      best_object_source
    )
  )
}

mif_best <- readRDS(
  best_object_source
)

best_object_file <- file.path(
  "results",
  "best_mif2_object.rds"
)

saveRDS(
  mif_best,
  best_object_file
)

cat(
  "Combined results saved to: ",
  combined_file,
  "\nBest fit saved to: ",
  best_fit_file,
  "\nBest MIF2 object saved to: ",
  best_object_file,
  "\n\nBest fit:\n",
  sep = ""
)

print(
  best_fit,
  row.names = FALSE
)
