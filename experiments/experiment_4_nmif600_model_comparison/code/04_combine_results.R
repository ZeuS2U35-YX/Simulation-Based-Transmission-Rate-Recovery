# ============================================================
# Validate and combine raw Experiment 4 task outputs
#
# Usage:
# Rscript code/04_combine_results.R \
#   <model> <raw_results_root> <combined_output> [task_spec]
#
# model: gamma or constant
# task_spec examples: 1:200   or   1,50,100,150,200
# ============================================================

options(stringsAsFactors = FALSE)
source("config/experiment_config.R")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3L) {
  stop("Usage: Rscript code/04_combine_results.R <model> <raw_results_root> <combined_output> [task_spec]")
}

model_arg <- args[[1]]
raw_root <- args[[2]]
combined_dir <- args[[3]]
task_spec <- if (length(args) >= 4L) args[[4]] else paste0("1:", experiment_config$n_tasks)

if (!(model_arg %in% c("gamma", "constant"))) {
  stop("model must be 'gamma' or 'constant'.")
}
model_name <- if (model_arg == "gamma") "gamma_noise" else "constant_B"
accepted_model_names <- if (model_arg == "gamma") c("gamma_noise", "gamma_transition") else "constant_B"
expected_starts <- if (model_arg == "gamma") 9L else 6L

parse_task_spec <- function(x) {
  if (grepl("^[0-9]+:[0-9]+$", x)) {
    ends <- as.integer(strsplit(x, ":", fixed = TRUE)[[1]])
    return(seq.int(ends[[1]], ends[[2]]))
  }
  out <- suppressWarnings(as.integer(strsplit(x, ",", fixed = TRUE)[[1]]))
  if (any(is.na(out))) stop("Invalid task_spec: ", x)
  sort(unique(out))
}
expected_tasks <- parse_task_spec(task_spec)

dir.create(combined_dir, recursive = TRUE, showWarnings = FALSE)

required_files <- c(
  "COMPLETE", "mif2_results.csv", "evaluation_logliks.csv",
  "best_fit_summary.csv", "B_path.csv", "run_config.csv"
)

mif_list <- list()
eval_list <- list()
best_list <- list()
path_list <- list()
config_list <- list()
trace_list <- list()
problems <- character(0)

read_csv <- function(path) read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)

for (task_id in expected_tasks) {
  task_dir <- file.path(raw_root, sprintf("task_%03d", task_id))
  missing <- required_files[!file.exists(file.path(task_dir, required_files))]
  if (length(missing) > 0L) {
    problems <- c(problems, paste0("task_", sprintf("%03d", task_id), " missing: ", paste(missing, collapse = ", ")))
    next
  }

  mif <- read_csv(file.path(task_dir, "mif2_results.csv"))
  eval <- read_csv(file.path(task_dir, "evaluation_logliks.csv"))
  best <- read_csv(file.path(task_dir, "best_fit_summary.csv"))
  path <- read_csv(file.path(task_dir, "B_path.csv"))
  run_config <- read_csv(file.path(task_dir, "run_config.csv"))

  if (nrow(mif) != expected_starts) problems <- c(problems, paste0("task_", sprintf("%03d", task_id), " has ", nrow(mif), " MIF rows; expected ", expected_starts))
  if (nrow(best) != 1L) problems <- c(problems, paste0("task_", sprintf("%03d", task_id), " best summary has ", nrow(best), " rows"))
  if (nrow(eval) != expected_starts * experiment_config$n_pf_evals) problems <- c(problems, paste0("task_", sprintf("%03d", task_id), " evaluation table has unexpected row count"))

  if (!all(mif$task_id == task_id) || !all(best$task_id == task_id) || (nrow(path) > 0L && !all(path$task_id == task_id))) {
    problems <- c(problems, paste0("task_", sprintf("%03d", task_id), " contains a task_id mismatch"))
  }
  if (!all(mif$model %in% accepted_model_names) ||
      !all(best$model %in% accepted_model_names) ||
      (nrow(path) > 0L && !all(path$model %in% accepted_model_names))) {
    problems <- c(problems, paste0("task_", sprintf("%03d", task_id), " contains a model mismatch"))
  }
  if (!all(mif$Nmif == experiment_config$Nmif) || !all(best$Nmif == experiment_config$Nmif)) {
    problems <- c(problems, paste0("task_", sprintf("%03d", task_id), " is not Nmif=", experiment_config$Nmif))
  }
  if (best$status[[1]] != "success") {
    problems <- c(problems, paste0("task_", sprintf("%03d", task_id), " best-fit status is ", best$status[[1]]))
  }
  if (nrow(path) != 70L) {
    problems <- c(problems, paste0("task_", sprintf("%03d", task_id), " has ", nrow(path), " B-path rows; expected 70"))
  }

  # Canonicalize the machine-readable label while accepting legacy raw
  # outputs produced before the Gamma-noise terminology update.
  mif$model <- model_name
  eval$model <- model_name
  best$model <- model_name
  if (nrow(path) > 0L) path$model <- model_name
  if ("setting" %in% names(run_config) && "value" %in% names(run_config)) {
    run_config$value[run_config$setting == "model"] <- model_name
  }

  mif_list[[as.character(task_id)]] <- mif
  eval_list[[as.character(task_id)]] <- eval
  best_list[[as.character(task_id)]] <- best
  path_list[[as.character(task_id)]] <- path
  run_config$task_id <- task_id
  run_config$model <- model_name
  config_list[[as.character(task_id)]] <- run_config

  trace_file <- file.path(task_dir, "mif2_traces.csv")
  if (file.exists(trace_file)) {
    trace <- read_csv(trace_file)
    trace$model <- model_name
    trace_list[[as.character(task_id)]] <- trace
  }
}

bind_or_empty <- function(x) {
  if (length(x) == 0L) data.frame() else do.call(rbind, x)
}

combined_mif <- bind_or_empty(mif_list)
combined_eval <- bind_or_empty(eval_list)
combined_best <- bind_or_empty(best_list)
combined_path <- bind_or_empty(path_list)
combined_config <- bind_or_empty(config_list)
combined_trace <- bind_or_empty(trace_list)

write.csv(combined_mif, file.path(combined_dir, "combined_mif2_results.csv"), row.names = FALSE)
write.csv(combined_eval, file.path(combined_dir, "combined_evaluation_logliks.csv"), row.names = FALSE)
write.csv(combined_best, file.path(combined_dir, "combined_best_fit_summary.csv"), row.names = FALSE)
write.csv(combined_path, file.path(combined_dir, "combined_B_paths.csv"), row.names = FALSE)
write.csv(combined_config, file.path(combined_dir, "combined_run_config.csv"), row.names = FALSE)
if (nrow(combined_trace) > 0L) {
  write.csv(combined_trace, file.path(combined_dir, "combined_mif2_traces.csv"), row.names = FALSE)
}

found_tasks <- if (nrow(combined_best) > 0L) sort(unique(combined_best$task_id)) else integer(0)
run_check <- data.frame(
  model = model_name,
  task_spec = task_spec,
  n_expected_tasks = length(expected_tasks),
  n_found_tasks = length(found_tasks),
  n_missing_tasks = length(setdiff(expected_tasks, found_tasks)),
  n_problem_messages = length(problems),
  all_task_ids_unique = nrow(combined_best) == length(unique(combined_best$task_id)),
  all_nmif_600 = nrow(combined_best) > 0L && all(combined_best$Nmif == experiment_config$Nmif),
  all_status_success = nrow(combined_best) > 0L && all(combined_best$status == "success"),
  stringsAsFactors = FALSE
)
write.csv(run_check, file.path(combined_dir, "run_check_summary.csv"), row.names = FALSE)
writeLines(problems, file.path(combined_dir, "combine_problems.txt"))

if (length(problems) > 0L) {
  stop(
    "Combination found ", length(problems),
    " problem(s). See ", file.path(combined_dir, "combine_problems.txt"), "."
  )
}

cat("Combined ", length(expected_tasks), " ", model_name, " tasks into ", combined_dir, ".\n", sep = "")
