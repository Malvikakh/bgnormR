#' Check if required packages are available
#'
#' Helper function to check if required packages are installed and loaded.
#'
#' @param packages Character vector of package names to check.
#' @param func_name Name of the calling function (for error messages).
#'
#' @return Invisible TRUE if all packages are available.
#'
#' @keywords internal
check_required_packages <- function(packages, func_name = NULL) {
    missing_pkgs <- character(0)
    
    for (pkg in packages) {
        if (!requireNamespace(pkg, quietly = TRUE)) {
            missing_pkgs <- c(missing_pkgs, pkg)
        }
    }
    
    if (length(missing_pkgs) > 0) {
        func_msg <- if (!is.null(func_name)) paste0(" for ", func_name) else ""
        stop("Required package(s) not installed", func_msg, ": ",
             paste(missing_pkgs, collapse = ", "),
             "\nInstall with: install.packages(c('", 
             paste(missing_pkgs, collapse = "', '"), "'))",
             call. = FALSE)
    }
    
    invisible(TRUE)
}


#' Validate file paths
#'
#' Helper function to validate that files exist and are readable.
#'
#' @param files Character vector of file paths to validate.
#' @param file_type Optional description of file type for error messages.
#'
#' @return Invisible TRUE if all files are valid.
#'
#' @keywords internal
validate_files <- function(files, file_type = "File") {
    for (f in files) {
        if (!file.exists(f)) {
            stop(file_type, " does not exist: ", f, call. = FALSE)
        }
        if (!file.access(f, mode = 4) == 0) {
            stop(file_type, " is not readable: ", f, call. = FALSE)
        }
    }
    invisible(TRUE)
}


#' Validate numeric parameters
#'
#' Helper function to validate numeric parameters are within acceptable ranges.
#'
#' @param value The value to validate.
#' @param name Parameter name (for error messages).
#' @param min Minimum acceptable value (inclusive). NULL for no minimum.
#' @param max Maximum acceptable value (inclusive). NULL for no maximum.
#' @param allow_negative Logical, whether negative values are allowed.
#'
#' @return Invisible TRUE if validation passes.
#'
#' @keywords internal
validate_numeric <- function(value, name, min = NULL, max = NULL, 
                            allow_negative = FALSE) {
    if (!is.numeric(value)) {
        stop(name, " must be numeric", call. = FALSE)
    }
    
    if (length(value) != 1) {
        stop(name, " must be a single value", call. = FALSE)
    }
    
    if (!allow_negative && value < 0) {
        stop(name, " must be non-negative", call. = FALSE)
    }
    
    if (!is.null(min) && value < min) {
        stop(name, " must be >= ", min, call. = FALSE)
    }
    
    if (!is.null(max) && value > max) {
        stop(name, " must be <= ", max, call. = FALSE)
    }
    
    invisible(TRUE)
}


#' Validate matrix or array dimensions
#'
#' Helper function to validate matrix/array dimensions match expected values.
#'
#' @param x Matrix or array to validate.
#' @param expected_dims Expected dimensions. NULL entries are not checked.
#' @param name Object name (for error messages).
#'
#' @return Invisible TRUE if validation passes.
#'
#' @keywords internal
validate_dimensions <- function(x, expected_dims = NULL, name = "Object") {
    if (!is.matrix(x) && !is.array(x)) {
        stop(name, " must be a matrix or array", call. = FALSE)
    }
    
    if (!is.null(expected_dims)) {
        actual_dims <- dim(x)
        if (length(expected_dims) != length(actual_dims)) {
            stop(name, " has wrong number of dimensions: expected ",
                 length(expected_dims), ", got ", length(actual_dims),
                 call. = FALSE)
        }
        
        for (i in seq_along(expected_dims)) {
            if (!is.null(expected_dims[i]) && expected_dims[i] != actual_dims[i]) {
                stop(name, " dimension ", i, " mismatch: expected ",
                     expected_dims[i], ", got ", actual_dims[i],
                     call. = FALSE)
            }
        }
    }
    
    invisible(TRUE)
}


#' Create a simple progress bar
#'
#' Helper function to create and update a text-based progress bar.
#'
#' @param current Current iteration number.
#' @param total Total number of iterations.
#' @param width Width of the progress bar in characters (default: 50).
#' @param prefix Optional prefix text.
#'
#' @return Invisible NULL. Prints progress bar to console.
#'
#' @keywords internal
show_progress <- function(current, total, width = 50, prefix = "") {
    percent <- current / total
    filled <- round(width * percent)
    bar <- paste0(
        "\r", prefix,
        "[", paste(rep("=", filled), collapse = ""),
        paste(rep(" ", width - filled), collapse = ""),
        "] ", round(percent * 100), "% (", current, "/", total, ")"
    )
    cat(bar)
    if (current == total) cat("\n")
    invisible(NULL)
}


#' Format bytes to human-readable size
#'
#' Helper function to convert bytes to human-readable format.
#'
#' @param bytes Number of bytes.
#'
#' @return Character string with formatted size.
#'
#' @keywords internal
format_bytes <- function(bytes) {
    units <- c("B", "KB", "MB", "GB", "TB")
    unit_idx <- 1
    size <- bytes
    
    while (size >= 1024 && unit_idx < length(units)) {
        size <- size / 1024
        unit_idx <- unit_idx + 1
    }
    
    paste0(round(size, 2), " ", units[unit_idx])
}


#' Get system information for diagnostics
#'
#' Gather system information useful for troubleshooting.
#'
#' @return Named list with system information.
#'
#' @examples
#' \dontrun{
#' info <- get_system_info()
#' print(info)
#' }
#'
#' @export
get_system_info <- function() {
    info <- list(
        r_version = R.version.string,
        platform = R.version$platform,
        os = Sys.info()["sysname"],
        memory_limit = if (.Platform$OS.type == "windows") memory.limit() else NA,
        java_version = tryCatch(
            system("java -version", intern = TRUE, ignore.stderr = FALSE),
            error = function(e) "Java not found or not configured"
        )
    )
    
    # Try to get Java heap size
    java_params <- getOption("java.parameters")
    if (!is.null(java_params)) {
        info$java_heap <- java_params
    } else {
        info$java_heap <- "Not set (default)"
    }
    
    class(info) <- c("bgnorm_sysinfo", "list")
    return(info)
}


#' Print system information
#'
#' @param x Object of class bgnorm_sysinfo.
#' @param ... Additional arguments (ignored).
#'
#' @export
print.bgnorm_sysinfo <- function(x, ...) {
    cat("=== BgNorm System Information ===\n")
    cat("R Version:     ", x$r_version, "\n")
    cat("Platform:      ", x$platform, "\n")
    cat("OS:            ", x$os, "\n")
    if (!is.na(x$memory_limit)) {
        cat("Memory Limit:  ", x$memory_limit, " MB\n")
    }
    cat("Java Heap:     ", x$java_heap, "\n")
    invisible(x)
}


#' Validate GMM summary data structure
#'
#' Check that GMM summary data has the required columns and structure.
#'
#' @param gmm_data Data frame to validate.
#'
#' @return Invisible TRUE if valid, stops with error otherwise.
#'
#' @keywords internal
validate_gmm_summary <- function(gmm_data) {
    required_cols <- c("sample", "channel", "component", "mean", "variance")
    
    if (!is.data.frame(gmm_data)) {
        stop("GMM summary data must be a data frame", call. = FALSE)
    }
    
    missing_cols <- setdiff(required_cols, colnames(gmm_data))
    if (length(missing_cols) > 0) {
        stop("GMM summary data is missing required columns: ",
             paste(missing_cols, collapse = ", "),
             call. = FALSE)
    }
    
    # Check for NA values in critical columns
    for (col in c("mean", "variance")) {
        if (any(is.na(gmm_data[[col]]))) {
            stop("GMM summary data contains NA values in column: ", col,
                 call. = FALSE)
        }
    }
    
    # Check variance is positive
    if (any(gmm_data$variance <= 0, na.rm = TRUE)) {
        stop("GMM summary data contains non-positive variance values",
             call. = FALSE)
    }
    
    invisible(TRUE)
}


#' Create summary statistics for GMM results
#'
#' Generate summary statistics from GMM processing results.
#'
#' @param gmm_result Result object from processChannelGMM.
#'
#' @return Data frame with summary statistics.
#'
#' @examples
#' \dontrun{
#' result <- processChannelGMM(...)
#' summary <- summarize_gmm_result(result)
#' print(summary)
#' }
#'
#' @export
summarize_gmm_result <- function(gmm_result) {
    if (!is.list(gmm_result)) {
        stop("gmm_result must be a list object from processChannelGMM",
             call. = FALSE)
    }
    
    summary_data <- data.frame(
        metric = character(),
        value = numeric(),
        stringsAsFactors = FALSE
    )
    
    # Add GMM summary info
    if ("gmm_summary" %in% names(gmm_result)) {
        gmm <- gmm_result$gmm_summary
        summary_data <- rbind(summary_data, data.frame(
            metric = "Number of components",
            value = max(gmm$component, na.rm = TRUE)
        ))
        
        summary_data <- rbind(summary_data, data.frame(
            metric = "Signal component",
            value = unique(gmm$G_sig)[1]
        ))
        
        summary_data <- rbind(summary_data, data.frame(
            metric = "Variance weight",
            value = unique(gmm$var_weight)[1]
        ))
        
        summary_data <- rbind(summary_data, data.frame(
            metric = "UQ75 normalization factor",
            value = unique(gmm$UQ75)[1]
        ))
    }
    
    # Add pixel data info
    if ("pixel_data" %in% names(gmm_result)) {
        pixels <- gmm_result$pixel_data
        summary_data <- rbind(summary_data, data.frame(
            metric = "Total pixels",
            value = nrow(pixels)
        ))
        
        summary_data <- rbind(summary_data, data.frame(
            metric = "Mean raw intensity",
            value = mean(pixels$raw_intensity, na.rm = TRUE)
        ))
        
        summary_data <- rbind(summary_data, data.frame(
            metric = "Mean adjusted intensity",
            value = mean(pixels$adj_intensity, na.rm = TRUE)
        ))
    }
    
    # Add cell data info
    if ("cell_data" %in% names(gmm_result) && !is.null(gmm_result$cell_data)) {
        cells <- gmm_result$cell_data
        summary_data <- rbind(summary_data, data.frame(
            metric = "Number of cells",
            value = nrow(cells)
        ))
        
        summary_data <- rbind(summary_data, data.frame(
            metric = "Mean cell median (adjusted)",
            value = mean(cells$median_adj, na.rm = TRUE)
        ))
    }
    
    return(summary_data)
}
