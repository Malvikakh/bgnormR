#' Complete GMM normalization workflow
#'
#' High-level wrapper function that performs the complete GMM normalization
#' workflow from raw image data to QC visualization.
#'
#' @param image_file Path to the multi-channel image file.
#' @param channel_names Character vector of channel names.
#' @param cell_mask Path to cell segmentation mask RDS file.
#' @param output_dir Directory for output files.
#' @param cofactor Cofactor for log transformation (default: 150).
#' @param bit_depth Bit depth of input image (default: 16).
#' @param G Number of GMM components (default: 3).
#' @param create_qc Logical, whether to create QC heatmaps (default: TRUE).
#' @param save_results Logical, whether to save results to files (default: TRUE).
#' @param verbose Logical, whether to print progress messages (default: TRUE).
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{gmm_results}: List of GMM results for each channel
#'   \item \code{gmm_summary}: Combined GMM summary data frame
#'   \item \code{cell_data}: Combined cell-level data
#'   \item \code{qc_plots}: Cohen's d QC heatmap plot (if create_qc = TRUE)
#' }
#'
#' @examples
#' \dontrun{
#' # Run complete workflow
#' results <- run_gmm_workflow(
#'   image_file = "sample.tif",
#'   channel_names = c("CD3", "CD8", "CD20"),
#'   cell_mask = "segmentation.rds",
#'   output_dir = "results"
#' )
#' 
#' # View QC heatmap
#' print(results$qc_plots)
#' 
#' # Access processed data
#' head(results$cell_data)
#' }
#'
#' @export
run_gmm_workflow <- function(image_file, channel_names, cell_mask, output_dir,
                             cofactor = 150, bit_depth = 16, G = 3,
                             create_qc = TRUE, save_results = TRUE,
                             verbose = TRUE) {
    
    # Validate inputs
    validate_files(c(image_file, cell_mask), "Input file")
    validate_numeric(cofactor, "cofactor", min = 1)
    validate_numeric(bit_depth, "bit_depth", min = 1, max = 32)
    validate_numeric(G, "G", min = 1, max = 10)
    
    if (!dir.exists(output_dir)) {
        dir.create(output_dir, recursive = TRUE)
        if (verbose) message("Created output directory: ", output_dir)
    }
    
    # Step 1: Process image channels
    if (verbose) message("Step 1/4: Processing image channels...")
    processed_dir <- file.path(output_dir, "processed")
    
    processMultiChannelImage(
        file = image_file,
        channel_names = channel_names,
        output_dir = processed_dir,
        cofactor = cofactor,
        bit_depth = bit_depth,
        parallel = FALSE
    )
    
    # Step 2: Apply GMM normalization to each channel
    if (verbose) message("Step 2/4: Applying GMM normalization...")
    gmm_results <- list()
    base_name <- tools::file_path_sans_ext(basename(image_file))
    
    for (i in seq_along(channel_names)) {
        ch <- channel_names[i]
        if (verbose) show_progress(i, length(channel_names), prefix = "  Channel: ")
        
        # Find the processed file
        processed_file <- list.files(
            processed_dir,
            pattern = paste0(ch, "_log2\\.rds$"),
            full.names = TRUE
        )[1]
        
        if (is.na(processed_file)) {
            warning("Could not find processed file for channel: ", ch)
            next
        }
        
        result <- processChannelGMM(
            file = processed_file,
            channel_name = ch,
            cell_mask = cell_mask,
            cofactor = cofactor,
            bit_depth = bit_depth,
            G = G
        )
        
        result$gmm_summary$sample <- base_name
        if (!is.null(result$cell_data)) {
            result$cell_data$sample <- base_name
        }
        
        gmm_results[[ch]] <- result
    }
    
    # Step 3: Combine results
    if (verbose) message("Step 3/4: Combining results...")
    gmm_summary <- do.call(rbind, lapply(gmm_results, function(x) x$gmm_summary))
    cell_data <- do.call(rbind, lapply(gmm_results, function(x) x$cell_data))
    
    # Step 4: Create QC visualizations
    qc_plots <- NULL
    if (create_qc && check_required_packages(c("ggplot2", "dplyr"), NULL)) {
        if (verbose) message("Step 4/4: Creating QC visualizations...")
        
        qc_plots <- createQCHeatmap(gmm_summary)
    } else if (verbose) {
        message("Step 4/4: Skipping QC visualizations (packages not available)")
    }
    
    # Save results
    if (save_results) {
        if (verbose) message("Saving results to disk...")
        saveRDS(gmm_results, file.path(output_dir, "gmm_results.rds"))
        saveRDS(gmm_summary, file.path(output_dir, "gmm_summary.rds"))
        write.csv(gmm_summary, file.path(output_dir, "gmm_summary.csv"), row.names = FALSE)
        
        if (!is.null(cell_data)) {
            saveRDS(cell_data, file.path(output_dir, "cell_data.rds"))
            write.csv(cell_data, file.path(output_dir, "cell_data.csv"), row.names = FALSE)
        }
        
        if (!is.null(qc_plots)) {
            if (requireNamespace("ggplot2", quietly = TRUE)) {
                ggplot2::ggsave(
                    file.path(output_dir, "qc_cohens_d.png"),
                    qc_plots, width = 10, height = 6
                )
            }
        }
    }
    
    if (verbose) message("Workflow complete!")
    
    # Return results
    return(structure(list(
        gmm_results = gmm_results,
        gmm_summary = gmm_summary,
        cell_data = cell_data,
        qc_plots = qc_plots
    ), class = c("bgnorm_workflow_result", "list")))
}


#' Print method for workflow results
#'
#' @param x Object of class bgnorm_workflow_result.
#' @param ... Additional arguments (ignored).
#'
#' @export
print.bgnorm_workflow_result <- function(x, ...) {
    cat("=== BgNorm Workflow Results ===\n")
    cat("Channels processed:  ", length(x$gmm_results), "\n")
    cat("Samples:             ", length(unique(x$gmm_summary$sample)), "\n")
    if (!is.null(x$cell_data)) {
        cat("Total cells:         ", nrow(x$cell_data), "\n")
    }
    cat("QC plots available:  ", !is.null(x$qc_plots), "\n")
    cat("\nUse summary() for detailed statistics\n")
    invisible(x)
}


#' Summary method for workflow results
#'
#' @param object Object of class bgnorm_workflow_result.
#' @param ... Additional arguments (ignored).
#'
#' @export
summary.bgnorm_workflow_result <- function(object, ...) {
    cat("=== BgNorm Workflow Summary ===\n\n")
    
    cat("Channels:\n")
    for (ch in names(object$gmm_results)) {
        cat("  -", ch, "\n")
        result_summary <- summarize_gmm_result(object$gmm_results[[ch]])
        for (i in 1:nrow(result_summary)) {
            cat("    ", result_summary$metric[i], ": ", 
                round(result_summary$value[i], 3), "\n", sep = "")
        }
    }
    
    cat("\nOverall Statistics:\n")
    cat("  Total GMM components fitted:", nrow(object$gmm_summary), "\n")
    if (!is.null(object$cell_data)) {
        cat("  Cells analyzed:", nrow(object$cell_data), "\n")
        cat("  Mean cell intensity (adjusted):", 
            round(mean(object$cell_data$median_adj, na.rm = TRUE), 3), "\n")
    }
    
    invisible(object)
}

