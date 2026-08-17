# ============================================================
# File, seed, and task-directory helpers for Experiment 4
# ============================================================

read_paramlist_task <- function(paramlist_file, task_id, required_columns) {
  if (!file.exists(paramlist_file)) {
    stop("Could not find paramlist file: ", paramlist_file)
  }

  paramlist <- read.csv(
    paramlist_file,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  missing_columns <- setdiff(required_columns, names(paramlist))
  if (length(missing_columns) > 0) {
    stop(
      "paramlist.csv is missing column(s): ",
      paste(missing_columns, collapse = ", ")
    )
  }

  task_row <- paramlist[paramlist$task_id == task_id, , drop = FALSE]
  if (nrow(task_row) != 1L) {
    stop(
      "Expected exactly one paramlist row for task_id = ", task_id,
      "; found ", nrow(task_row), "."
    )
  }

  task_row
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
      ". Remove or rename it before rerunning this task."
    )
  }

  job_id <- Sys.getenv("SLURM_JOB_ID", unset = "local")
  temp_dir <- file.path(
    output_root,
    sprintf(".task_%03d_tmp_%s_%s", task_id, job_id, Sys.getpid())
  )

  if (dir.exists(temp_dir)) {
    unlink(temp_dir, recursive = TRUE, force = TRUE)
  }
  dir.create(temp_dir, recursive = TRUE, showWarnings = FALSE)

  list(skip = FALSE, final_dir = final_dir, temp_dir = temp_dir)
}

commit_atomic_task <- function(temp_dir, final_dir) {
  writeLines(
    c(
      paste0("completed_at=", format(Sys.time(), tz = "UTC", usetz = TRUE)),
      paste0("host=", Sys.info()[["nodename"]]),
      paste0("slurm_job_id=", Sys.getenv("SLURM_JOB_ID", unset = "local")),
      paste0("slurm_array_task_id=", Sys.getenv("SLURM_ARRAY_TASK_ID", unset = "local"))
    ),
    file.path(temp_dir, "COMPLETE")
  )

  if (!file.rename(temp_dir, final_dir)) {
    stop("Could not atomically rename ", temp_dir, " to ", final_dir, ".")
  }
}

file_md5 <- function(path) {
  if (!file.exists(path)) {
    stop("Cannot calculate MD5; file does not exist: ", path)
  }
  unname(tools::md5sum(path)[[1]])
}

read_shared_data <- function(shared_data_root, task_id) {
  task_dir <- file.path(shared_data_root, sprintf("task_%03d", task_id))
  required_files <- c(
    "COMPLETE",
    "observed_data.csv",
    "simulated_data.csv",
    "simulation_metadata.csv",
    "data_checksums.csv"
  )

  missing <- required_files[!file.exists(file.path(task_dir, required_files))]
  if (length(missing) > 0) {
    stop(
      "Shared data for task ", task_id, " are incomplete. Missing: ",
      paste(missing, collapse = ", ")
    )
  }

  checksums <- read.csv(
    file.path(task_dir, "data_checksums.csv"),
    stringsAsFactors = FALSE
  )

  observed_file <- file.path(task_dir, "observed_data.csv")
  simulated_file <- file.path(task_dir, "simulated_data.csv")

  expected_observed_md5 <- checksums$md5[checksums$file == "observed_data.csv"]
  expected_simulated_md5 <- checksums$md5[checksums$file == "simulated_data.csv"]

  if (length(expected_observed_md5) != 1L || length(expected_simulated_md5) != 1L) {
    stop("Invalid data_checksums.csv for task ", task_id, ".")
  }

  actual_observed_md5 <- file_md5(observed_file)
  actual_simulated_md5 <- file_md5(simulated_file)

  if (!identical(actual_observed_md5, expected_observed_md5)) {
    stop("Observed-data MD5 mismatch for task ", task_id, ".")
  }
  if (!identical(actual_simulated_md5, expected_simulated_md5)) {
    stop("Simulated-data MD5 mismatch for task ", task_id, ".")
  }

  list(
    task_dir = task_dir,
    observed_data = read.csv(observed_file, check.names = FALSE),
    simulated_data = read.csv(simulated_file, check.names = FALSE),
    metadata = read.csv(
      file.path(task_dir, "simulation_metadata.csv"),
      check.names = FALSE,
      stringsAsFactors = FALSE
    ),
    observed_data_md5 = actual_observed_md5,
    simulated_data_md5 = actual_simulated_md5
  )
}

safe_trace_data <- function(mif_object, task_id, model, run, start_values) {
  out <- tryCatch(
    as.data.frame(traces(mif_object)),
    error = function(e) NULL
  )

  if (is.null(out)) {
    return(NULL)
  }

  out$mif_iteration <- seq_len(nrow(out)) - 1L
  out$task_id <- task_id
  out$model <- model
  out$run <- run

  for (name in names(start_values)) {
    out[[paste0("start_", name)]] <- start_values[[name]]
  }

  out
}

write_run_config <- function(path, values) {
  values <- c(
    values,
    list(
      R_version = R.version.string,
      pomp_version = as.character(packageVersion("pomp")),
      hostname = Sys.info()[["nodename"]],
      slurm_job_id = Sys.getenv("SLURM_JOB_ID", unset = "local"),
      slurm_array_task_id = Sys.getenv("SLURM_ARRAY_TASK_ID", unset = "local"),
      written_at_utc = format(Sys.time(), tz = "UTC", usetz = TRUE)
    )
  )

  data <- data.frame(
    setting = names(values),
    value = vapply(values, function(x) paste(x, collapse = ","), character(1)),
    stringsAsFactors = FALSE
  )

  write.csv(data, path, row.names = FALSE)
}
