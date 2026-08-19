# ============================================================
# I/O and validation helpers for the paired Experiment 5 batch
# ============================================================

file_md5 <- function(path) {
  if (!file.exists(path)) {
    stop("Cannot calculate MD5; file does not exist: ", path)
  }
  unname(tools::md5sum(path)[[1]])
}

read_csv <- function(path) {
  if (!file.exists(path)) stop("Required file does not exist: ", path)
  read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
}

require_exact_task_ids <- function(x, n_tasks, label) {
  expected <- seq_len(n_tasks)
  actual <- sort(as.integer(x))
  if (length(actual) != n_tasks || anyDuplicated(actual) ||
      !identical(actual, expected)) {
    stop(
      label, " must contain each task_id from 1 through ", n_tasks,
      " exactly once."
    )
  }
  invisible(TRUE)
}

read_exp4_shared_task <- function(shared_data_root, task_id) {
  task_dir <- file.path(shared_data_root, sprintf("task_%03d", task_id))
  required_files <- c(
    "COMPLETE", "observed_data.csv", "simulated_data.csv",
    "simulation_metadata.csv", "data_checksums.csv"
  )
  missing <- required_files[!file.exists(file.path(task_dir, required_files))]
  if (length(missing) > 0L) {
    stop(
      "Experiment 4 shared data for task ", task_id,
      " are incomplete. Missing: ", paste(missing, collapse = ", ")
    )
  }

  checksums <- read_csv(file.path(task_dir, "data_checksums.csv"))
  observed_path <- file.path(task_dir, "observed_data.csv")
  simulated_path <- file.path(task_dir, "simulated_data.csv")
  expected_observed_md5 <- checksums$md5[
    checksums$file == "observed_data.csv"
  ]
  expected_simulated_md5 <- checksums$md5[
    checksums$file == "simulated_data.csv"
  ]
  if (length(expected_observed_md5) != 1L ||
      length(expected_simulated_md5) != 1L) {
    stop("Invalid Experiment 4 checksum table for task ", task_id, ".")
  }

  observed_md5 <- file_md5(observed_path)
  if (!identical(observed_md5, expected_observed_md5)) {
    stop("Observed-data checksum mismatch for task ", task_id, ".")
  }

  simulated_md5 <- file_md5(simulated_path)
  if (!identical(simulated_md5, expected_simulated_md5)) {
    stop("Simulated-data checksum mismatch for task ", task_id, ".")
  }

  metadata <- read_csv(file.path(task_dir, "simulation_metadata.csv"))
  required_metadata_columns <- c(
    "task_id", "accepted", "acceptance_threshold", "max_H"
  )
  missing_metadata_columns <- setdiff(
    required_metadata_columns, names(metadata)
  )
  if (nrow(metadata) != 1L || length(missing_metadata_columns) > 0L) {
    stop(
      "Invalid Experiment 4 simulation metadata for task ", task_id,
      if (length(missing_metadata_columns) > 0L) {
        paste0(
          ". Missing: ", paste(missing_metadata_columns, collapse = ", ")
        )
      } else {
        "."
      }
    )
  }
  if (as.integer(metadata$task_id[[1]]) != task_id) {
    stop("Simulation-metadata task_id mismatch for task ", task_id, ".")
  }

  accepted <- as.logical(metadata$accepted[[1]])
  if (length(accepted) != 1L || is.na(accepted) || !identical(accepted, TRUE)) {
    stop("Task ", task_id, " does not have accepted=TRUE.")
  }

  acceptance_threshold <- as.numeric(metadata$acceptance_threshold[[1]])
  if (length(acceptance_threshold) != 1L ||
      !is.finite(acceptance_threshold) ||
      !identical(acceptance_threshold, 20)) {
    stop("Task ", task_id, " does not have acceptance_threshold=20.")
  }

  simulated_data <- read_csv(simulated_path)
  if (!"H" %in% names(simulated_data)) {
    stop("Simulated data for task ", task_id, " do not contain column H.")
  }
  if (!is.numeric(simulated_data$H) || length(simulated_data$H) == 0L ||
      !all(is.finite(simulated_data$H))) {
    stop("Simulated-data H values for task ", task_id, " are not all finite.")
  }
  recomputed_max_H <- max(simulated_data$H)
  if (!(recomputed_max_H > acceptance_threshold)) {
    stop(
      "Task ", task_id,
      " fails the independently recomputed max_H acceptance condition."
    )
  }

  metadata_max_H <- as.numeric(metadata$max_H[[1]])
  if (length(metadata_max_H) != 1L || !is.finite(metadata_max_H) ||
      !identical(metadata_max_H, as.numeric(recomputed_max_H))) {
    stop(
      "Simulation-metadata max_H does not exactly match recomputed max_H ",
      "for task ", task_id, "."
    )
  }

  observed_data <- read_csv(observed_path)
  if (nrow(observed_data) != 70L ||
      !all(c("week", "reports") %in% names(observed_data)) ||
      is.unsorted(observed_data$week, strictly = TRUE) ||
      anyDuplicated(observed_data$week)) {
    stop("Observed data for task ", task_id, " are not the expected 70 rows.")
  }

  list(
    task_dir = task_dir,
    observed_data = observed_data,
    simulated_data = simulated_data,
    metadata = metadata,
    observed_data_md5 = observed_md5,
    simulated_data_md5 = simulated_md5,
    acceptance_threshold = acceptance_threshold,
    recomputed_max_H = recomputed_max_H
  )
}

read_manifest_task <- function(manifest_file, task_id, n_tasks) {
  manifest <- read_csv(manifest_file)
  required_columns <- c(
    "task_id", "simulation_seed", "simulation_attempt",
    "observed_data_md5", "simulated_data_md5", "acceptance_threshold",
    "recomputed_max_H", "accepted"
  )
  missing_columns <- setdiff(required_columns, names(manifest))
  if (length(missing_columns) > 0L) {
    stop(
      "Input manifest is missing: ", paste(missing_columns, collapse = ", ")
    )
  }
  require_exact_task_ids(manifest$task_id, n_tasks, "Input manifest")
  row <- manifest[manifest$task_id == task_id, , drop = FALSE]
  if (nrow(row) != 1L) {
    stop("Input manifest does not have exactly one row for task ", task_id, ".")
  }
  row
}

start_atomic_task <- function(output_root, task_id) {
  dir.create(output_root, recursive = TRUE, showWarnings = FALSE)
  final_dir <- file.path(output_root, sprintf("task_%03d", task_id))
  complete_file <- file.path(final_dir, "COMPLETE")

  if (file.exists(complete_file)) {
    cat("Task output is already complete; skipping: ", final_dir, "\n", sep = "")
    return(list(skip = TRUE, final_dir = final_dir, temp_dir = NA_character_))
  }
  if (dir.exists(final_dir)) {
    stop(
      "An incomplete task directory already exists: ", final_dir,
      ". Move it aside before rerunning this task."
    )
  }

  job_id <- Sys.getenv("SLURM_JOB_ID", unset = "local")
  temp_dir <- file.path(
    output_root,
    sprintf(".task_%03d_tmp_%s_%s", task_id, job_id, Sys.getpid())
  )
  if (dir.exists(temp_dir)) unlink(temp_dir, recursive = TRUE, force = TRUE)
  dir.create(temp_dir, recursive = TRUE, showWarnings = FALSE)

  list(skip = FALSE, final_dir = final_dir, temp_dir = temp_dir)
}

commit_atomic_task <- function(temp_dir, final_dir) {
  writeLines(
    c(
      paste0("completed_at=", format(Sys.time(), tz = "UTC", usetz = TRUE)),
      paste0("host=", Sys.info()[["nodename"]]),
      paste0("slurm_job_id=", Sys.getenv("SLURM_JOB_ID", unset = "local")),
      paste0(
        "slurm_array_task_id=",
        Sys.getenv("SLURM_ARRAY_TASK_ID", unset = "local")
      )
    ),
    file.path(temp_dir, "COMPLETE")
  )
  if (!file.rename(temp_dir, final_dir)) {
    stop("Could not atomically rename ", temp_dir, " to ", final_dir, ".")
  }
}

write_run_config <- function(path, values) {
  values <- c(
    values,
    list(
      R_version = R.version.string,
      pomp_version = as.character(packageVersion("pomp")),
      hostname = Sys.info()[["nodename"]],
      slurm_job_id = Sys.getenv("SLURM_JOB_ID", unset = "local"),
      slurm_array_task_id = Sys.getenv(
        "SLURM_ARRAY_TASK_ID", unset = "local"
      ),
      written_at_utc = format(Sys.time(), tz = "UTC", usetz = TRUE)
    )
  )
  out <- data.frame(
    setting = names(values),
    value = vapply(values, function(x) paste(x, collapse = ","), character(1)),
    stringsAsFactors = FALSE
  )
  write.csv(out, path, row.names = FALSE)
}

safe_trace_data <- function(mif_object, task_id, run, start_values) {
  out <- tryCatch(
    as.data.frame(pomp::traces(mif_object)),
    error = function(e) NULL
  )
  if (is.null(out)) return(NULL)

  out$mif_iteration <- seq_len(nrow(out)) - 1L
  out$task_id <- task_id
  out$model <- "bspline_B"
  out$run <- run
  for (name in names(start_values)) {
    out[[paste0("start_", name)]] <- start_values[[name]]
  }
  out
}
