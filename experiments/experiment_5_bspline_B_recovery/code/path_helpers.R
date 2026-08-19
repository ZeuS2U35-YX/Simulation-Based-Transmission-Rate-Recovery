get_script_directory <- function() {
  command_args <- commandArgs(trailingOnly = FALSE)
  file_argument <- grep("^--file=", command_args, value = TRUE)

  if (length(file_argument) > 0L) {
    return(dirname(normalizePath(
      sub("^--file=", "", file_argument[[1]]),
      mustWork = FALSE
    )))
  }

  source_file <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)

  if (!is.null(source_file)) {
    return(dirname(normalizePath(source_file, mustWork = FALSE)))
  }

  normalizePath(getwd(), mustWork = FALSE)
}

get_experiment_directory <- function() {
  script_directory <- get_script_directory()

  if (basename(script_directory) == "code") {
    dirname(script_directory)
  } else {
    script_directory
  }
}

ensure_experiment_directories <- function(experiment_directory) {
  directories <- file.path(
    experiment_directory,
    c("data", "results", "figures")
  )

  invisible(lapply(
    directories,
    dir.create,
    recursive = TRUE,
    showWarnings = FALSE
  ))
}
