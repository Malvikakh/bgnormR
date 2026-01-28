#' Process a single channel from a multi-channel image
#'
#' Reads and processes a single channel from a multi-channel TIFF image file,
#' applying median filtering and log transformation for background normalization.
#'
#' @param file Path to the input image file (must be readable by RBioFormats).
#' @param channel_index Numeric index of the channel to process (1-based indexing).
#' @param channel_name Character string with the name/marker for this channel.
#' @param output_dir Directory where processed RDS files will be saved.
#' @param cofactor Numeric cofactor for log transformation (default: 16).
#' @param filter_size Size of the median filter kernel (default: 1).
#' @param bit_depth Bit depth of the input image (default: 8 for 8-bit images).
#' @param series Image series to read from multi-series files (default: 1).
#' @param resolution Image resolution level to read (default: 1).
#'
#' @return Invisibly returns the processed matrix. The function also saves the
#'   processed data as an RDS file.
#'
#' @details
#' The function performs the following steps:
#' \itemize{
#'   \item Reads the specified channel from the image file
#'   \item Applies median filtering using EBImage
#'   \item Performs log2 transformation: log2(x / cofactor + 1)
#'   \item Saves the result to an RDS file
#' }
#'
#' @examples
#' \dontrun{
#' # Process channel 1 from an image file
#' processChannel(
#'   file = "image.tif",
#'   channel_index = 1,
#'   channel_name = "DAPI",
#'   output_dir = "output"
#' )
#' }
#'
#' @importFrom RBioFormats read.image
#' @importFrom EBImage medianFilter
#' @export
processChannel <- function(file, channel_index, channel_name, output_dir,
                          cofactor = 16, filter_size = 1, bit_depth = 8,
                          series = 1, resolution = 1) {
    
    # Validate inputs
    if (!file.exists(file)) {
        stop("Input file does not exist: ", file)
    }
    
    if (!dir.exists(output_dir)) {
        dir.create(output_dir, recursive = TRUE)
    }
    
    # Extract base filename
    base_filename <- sub(".* ", "", tools::file_path_sans_ext(basename(file)))
    
    message(paste("Processing channel", channel_name, "(index:", channel_index, ")"))
    
    # Read the image for the current channel
    img <- RBioFormats::read.image(
        file,
        subset = list(C = channel_index),
        series = series,
        resolution = resolution,
        normalize = FALSE
    )
    
    x <- img@.Data
    
    # Check pixel value range
    pixel_range <- range(x, na.rm = TRUE)
    message(paste("Pixel value range:", pixel_range[1], "to", pixel_range[2]))
    
    # Convert to matrix and apply median filter
    max_value <- 2^bit_depth - 1
    x <- as.matrix(EBImage::medianFilter(x / max_value, filter_size)) * max_value
    
    # Log transformation
    x <- log2(x / cofactor + 1)
    
    # Save processed data
    output_file <- file.path(
        output_dir,
        paste0(base_filename, "_", channel_name, "_log2.rds")
    )
    saveRDS(x, file = output_file)
    message(paste("Saved to:", output_file))
    
    invisible(x)
}


#' Process all channels from a multi-channel image
#'
#' Reads metadata and processes all channels from a multi-channel TIFF image file.
#'
#' @param file Path to the input image file.
#' @param channel_names Character vector of channel/marker names. If NULL,
#'   channels will be named as "Channel_1", "Channel_2", etc.
#' @param output_dir Directory where processed RDS files will be saved.
#' @param cofactor Numeric cofactor for log transformation (default: 16).
#' @param filter_size Size of the median filter kernel (default: 1).
#' @param bit_depth Bit depth of the input image (default: 8).
#' @param series Image series to read (default: 1).
#' @param parallel Logical indicating whether to process channels in parallel
#'   (default: FALSE).
#' @param ncores Number of cores to use for parallel processing. If NULL,
#'   uses parallel::detectCores() - 1.
#'
#' @return A list of processed matrices, one per channel.
#'
#' @examples
#' \dontrun{
#' # Process all channels from an image
#' channels <- c("DAPI", "CD3", "CD8", "CD20")
#' results <- processMultiChannelImage(
#'   file = "image.tif",
#'   channel_names = channels,
#'   output_dir = "output"
#' )
#' }
#'
#' @importFrom RBioFormats read.metadata coreMetadata
#' @export
processMultiChannelImage <- function(file, channel_names = NULL, output_dir,
                                     cofactor = 16, filter_size = 1,
                                     bit_depth = 8, series = 1,
                                     parallel = FALSE, ncores = NULL) {
    
    # Read metadata to get number of channels
    metadata <- RBioFormats::read.metadata(file)
    num_channels <- RBioFormats::coreMetadata(metadata, series = series)$sizeC
    
    message(paste("Number of channels:", num_channels))
    
    # Generate channel names if not provided
    if (is.null(channel_names)) {
        channel_names <- paste0("Channel_", seq_len(num_channels))
    } else if (length(channel_names) != num_channels) {
        warning("Number of channel names does not match number of channels. ",
                "Using generic names.")
        channel_names <- paste0("Channel_", seq_len(num_channels))
    }
    
    # Function to process a single channel
    process_one <- function(j) {
        processChannel(
            file = file,
            channel_index = j,
            channel_name = channel_names[j],
            output_dir = output_dir,
            cofactor = cofactor,
            filter_size = filter_size,
            bit_depth = bit_depth,
            series = series
        )
    }
    
    # Process channels
    if (parallel) {
        if (!requireNamespace("parallel", quietly = TRUE)) {
            stop("Package 'parallel' is required for parallel processing")
        }
        
        if (is.null(ncores)) {
            ncores <- parallel::detectCores() - 1
        }
        
        message(paste("Processing channels in parallel using", ncores, "cores"))
        results <- parallel::mclapply(
            seq_len(num_channels),
            process_one,
            mc.cores = ncores
        )
    } else {
        message("Processing channels sequentially")
        results <- lapply(seq_len(num_channels), process_one)
    }
    
    names(results) <- channel_names
    invisible(results)
}


#' Batch process multiple image files
#'
#' Process multiple multi-channel image files in batch mode.
#'
#' @param file_pattern Pattern to match image files (e.g., "*.tif").
#' @param input_dir Directory containing input image files.
#' @param channel_names Character vector of channel/marker names, or path to
#'   a CSV file containing channel names (one per line).
#' @param output_dir Directory where processed RDS files will be saved.
#' @param cofactor Numeric cofactor for log transformation (default: 16).
#' @param filter_size Size of the median filter kernel (default: 1).
#' @param bit_depth Bit depth of input images (default: 8).
#' @param file_indices Optional numeric vector of file indices to process.
#'   If NULL, processes all files.
#' @param parallel Logical indicating whether to process channels in parallel
#'   within each image (default: FALSE).
#' @param ncores Number of cores for parallel processing within images.
#'
#' @return Invisibly returns a list of file paths that were processed.
#'
#' @examples
#' \dontrun{
#' # Process all TIFF files in a directory
#' batchProcessImages(
#'   file_pattern = "tif$",
#'   input_dir = "input_data",
#'   channel_names = "markers.csv",
#'   output_dir = "output_data"
#' )
#' }
#'
#' @export
batchProcessImages <- function(file_pattern, input_dir, channel_names,
                              output_dir, cofactor = 16, filter_size = 1,
                              bit_depth = 8, file_indices = NULL,
                              parallel = FALSE, ncores = NULL) {
    
    # Get list of image files
    image_files <- list.files(
        input_dir,
        pattern = file_pattern,
        full.names = TRUE
    )
    
    if (length(image_files) == 0) {
        stop("No files found matching pattern '", file_pattern,
             "' in directory '", input_dir, "'")
    }
    
    message(paste("Found", length(image_files), "image files"))
    
    # Read channel names if from file
    if (is.character(channel_names) && length(channel_names) == 1 &&
        file.exists(channel_names)) {
        message(paste("Reading channel names from:", channel_names))
        channel_names <- readLines(channel_names, warn = FALSE)
    }
    
    # Subset files if indices provided
    if (!is.null(file_indices)) {
        if (any(file_indices > length(image_files))) {
            stop("Some file indices are out of range")
        }
        image_files <- image_files[file_indices]
        message(paste("Processing", length(image_files), "selected files"))
    }
    
    # Process each file
    processed_files <- character(length(image_files))
    for (i in seq_along(image_files)) {
        file <- image_files[i]
        message(paste("\n=== Processing file", i, "of", length(image_files), "==="))
        message(paste("File:", basename(file)))
        
        processMultiChannelImage(
            file = file,
            channel_names = channel_names,
            output_dir = output_dir,
            cofactor = cofactor,
            filter_size = filter_size,
            bit_depth = bit_depth,
            parallel = parallel,
            ncores = ncores
        )
        
        processed_files[i] <- file
    }
    
    message(paste("\n=== Batch processing complete ==="))
    message(paste("Processed", length(processed_files), "files"))
    message(paste("Output directory:", output_dir))
    
    invisible(processed_files)
}
