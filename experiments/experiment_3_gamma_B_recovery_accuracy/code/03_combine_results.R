# ============================================================
# Combine results from the 200-task MIF2 simulation study
#
# Input:
#   Results/task_001/
#   ...
#   Results/task_200/
#
# Each task folder must contain:
#   - mif2_results.csv
#   - best_fit_summary.csv
#   - filtered_B_path.csv
#   - simulated_data.csv
#
# Output:
#   results/combined/
#     - combined_mif2_results.csv
#     - combined_best_fit_summary.csv
#     - combined_filtered_B_paths.csv
#     - combined_simulated_data.csv
#
# Usage:
# Rscript code/03_combine_results.R
#
# Optional:
# Rscript code/03_combine_results.R Results another_output_folder
# ============================================================

options(
  stringsAsFactors = FALSE
)

args <- commandArgs(
  trailingOnly = TRUE
)

if (length(args) >= 1) {
  results_folder <- args[[1]]
} else {
  results_folder <- "Results"
}

if (length(args) >= 2) {
  combined_folder <- args[[2]]
} else {
  combined_folder <- file.path("results", "combined")
}

if (!dir.exists(results_folder)) {
  stop(
    "Could not find results folder: ",
    results_folder
  )
}

dir.create(
  combined_folder,
  recursive = TRUE,
  showWarnings = FALSE
)

expected_task_ids <- seq_len(200L)

task_folders <- list.dirs(
  results_folder,
  recursive = FALSE,
  full.names = TRUE
)

task_folders <- task_folders[
  grepl(
    "^task_[0-9]+$",
    basename(task_folders)
  )
]

if (length(task_folders) == 0) {
  stop(
    "No task folders of the form task_001, task_002, ... ",
    "were found inside ",
    results_folder,
    "."
  )
}

task_ids_from_folders <- suppressWarnings(
  as.integer(
    sub(
      "^task_",
      "",
      basename(task_folders)
    )
  )
)

folder_order <- order(
  task_ids_from_folders
)

task_folders <- task_folders[
  folder_order
]

task_ids_from_folders <- task_ids_from_folders[
  folder_order
]

missing_task_ids <- setdiff(
  expected_task_ids,
  task_ids_from_folders
)

unexpected_task_ids <- setdiff(
  task_ids_from_folders,
  expected_task_ids
)

if (length(missing_task_ids) > 0) {
  stop(
    "Missing task folder(s): ",
    paste(
      sprintf(
        "task_%03d",
        missing_task_ids
      ),
      collapse = ", "
    )
  )
}

if (length(unexpected_task_ids) > 0) {
  warning(
    "Ignoring unexpected task folder(s): ",
    paste(
      sprintf(
        "task_%03d",
        unexpected_task_ids
      ),
      collapse = ", "
    )
  )
}

read_csv_checked <- function(
  file_path,
  description
) {

  if (!file.exists(file_path)) {
    stop(
      "Missing ",
      description,
      ": ",
      file_path
    )
  }

  read.csv(
    file_path,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}

mif_results_list <- vector(
  mode = "list",
  length = length(task_folders)
)

best_fit_list <- vector(
  mode = "list",
  length = length(task_folders)
)

B_path_list <- vector(
  mode = "list",
  length = length(task_folders)
)

simulated_data_list <- vector(
  mode = "list",
  length = length(task_folders)
)

for (i in seq_along(task_folders)) {

  task_folder <- task_folders[[i]]

  expected_task_id <- task_ids_from_folders[[i]]

  cat(
    "Reading ",
    basename(task_folder),
    "\n",
    sep = ""
  )

  mif_results_i <- read_csv_checked(
    file.path(
      task_folder,
      "mif2_results.csv"
    ),
    "mif2_results.csv"
  )

  best_fit_i <- read_csv_checked(
    file.path(
      task_folder,
      "best_fit_summary.csv"
    ),
    "best_fit_summary.csv"
  )

  B_path_i <- read_csv_checked(
    file.path(
      task_folder,
      "filtered_B_path.csv"
    ),
    "filtered_B_path.csv"
  )

  simulated_data_i <- read_csv_checked(
    file.path(
      task_folder,
      "simulated_data.csv"
    ),
    "simulated_data.csv"
  )

  if (nrow(mif_results_i) != 9L) {
    stop(
      basename(task_folder),
      " has ",
      nrow(mif_results_i),
      " mif2 rows, but 9 were expected."
    )
  }

  if (nrow(best_fit_i) != 1L) {
    stop(
      basename(task_folder),
      " has ",
      nrow(best_fit_i),
      " rows in best_fit_summary.csv, but 1 was expected."
    )
  }

  required_best_columns <- c(
    "task_id",
    "simulation_seed",
    "B0_hat",
    "sigma_beta_hat",
    "logLik",
    "logLik_se",
    "fit_success",
    "final_pf_success",
    "status"
  )

  missing_best_columns <- setdiff(
    required_best_columns,
    names(best_fit_i)
  )

  if (length(missing_best_columns) > 0) {
    stop(
      basename(task_folder),
      " is missing column(s) in best_fit_summary.csv: ",
      paste(
        missing_best_columns,
        collapse = ", "
      )
    )
  }

  saved_task_id <- as.integer(
    best_fit_i$task_id[[1]]
  )

  if (
    is.na(saved_task_id) ||
      saved_task_id != expected_task_id
  ) {
    stop(
      "Task ID mismatch in ",
      basename(task_folder),
      ". Folder suggests task ",
      expected_task_id,
      ", but best_fit_summary.csv says task ",
      saved_task_id,
      "."
    )
  }

  simulation_seed_i <- best_fit_i$simulation_seed[[1]]

  if (
    !("task_id" %in% names(simulated_data_i))
  ) {
    simulated_data_i$task_id <- expected_task_id
  }

  if (
    !("simulation_seed" %in% names(simulated_data_i))
  ) {
    simulated_data_i$simulation_seed <- simulation_seed_i
  }

  simulated_data_i <- simulated_data_i[
    ,
    c(
      "task_id",
      "simulation_seed",
      setdiff(
        names(simulated_data_i),
        c(
          "task_id",
          "simulation_seed"
        )
      )
    ),
    drop = FALSE
  ]

  mif_results_list[[i]] <- mif_results_i

  best_fit_list[[i]] <- best_fit_i

  B_path_list[[i]] <- B_path_i

  simulated_data_list[[i]] <- simulated_data_i
}

combined_mif_results <- do.call(
  rbind,
  mif_results_list
)

combined_best_fit_summary <- do.call(
  rbind,
  best_fit_list
)

combined_B_paths <- do.call(
  rbind,
  B_path_list
)

combined_simulated_data <- do.call(
  rbind,
  simulated_data_list
)

combined_mif_results <- combined_mif_results[
  order(
    combined_mif_results$task_id,
    combined_mif_results$run
  ),
  ,
  drop = FALSE
]

combined_best_fit_summary <- combined_best_fit_summary[
  order(
    combined_best_fit_summary$task_id
  ),
  ,
  drop = FALSE
]

combined_B_paths <- combined_B_paths[
  order(
    combined_B_paths$task_id,
    combined_B_paths$week
  ),
  ,
  drop = FALSE
]

combined_simulated_data <- combined_simulated_data[
  order(
    combined_simulated_data$task_id,
    combined_simulated_data$week
  ),
  ,
  drop = FALSE
]

write.csv(
  combined_mif_results,
  file.path(
    combined_folder,
    "combined_mif2_results.csv"
  ),
  row.names = FALSE
)

write.csv(
  combined_best_fit_summary,
  file.path(
    combined_folder,
    "combined_best_fit_summary.csv"
  ),
  row.names = FALSE
)

write.csv(
  combined_B_paths,
  file.path(
    combined_folder,
    "combined_filtered_B_paths.csv"
  ),
  row.names = FALSE
)

write.csv(
  combined_simulated_data,
  file.path(
    combined_folder,
    "combined_simulated_data.csv"
  ),
  row.names = FALSE
)

cat(
  "\nFinished combining results.\n",
  sep = ""
)

cat(
  "Number of tasks combined: ",
  length(task_folders),
  "\n",
  sep = ""
)

cat(
  "Rows in combined_mif2_results.csv: ",
  nrow(combined_mif_results),
  "\n",
  sep = ""
)

cat(
  "Rows in combined_best_fit_summary.csv: ",
  nrow(combined_best_fit_summary),
  "\n",
  sep = ""
)

cat(
  "Rows in combined_filtered_B_paths.csv: ",
  nrow(combined_B_paths),
  "\n",
  sep = ""
)

cat(
  "Rows in combined_simulated_data.csv: ",
  nrow(combined_simulated_data),
  "\n",
  sep = ""
)

cat(
  "\nSaved combined files in: ",
  combined_folder,
  "\n",
  sep = ""
)
