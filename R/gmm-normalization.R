#' Fit Gaussian Mixture Model to intensity data
#'
#' Fits a G-component Gaussian mixture model to log-transformed intensity values.
#'
#' @param intensity Numeric vector of intensity values (typically log-transformed).
#' @param G Number of mixture components (default: 3).
#'
#' @return A Mclust object containing the fitted GMM.
#'
#' @details
#' This function fits a Gaussian mixture model using the mclust package.
#' If fitting fails with G components, it automatically retries with G-1 components.
#'
#' @importFrom mclust Mclust
#' @export
fitGMMModel <- function(intensity, G = 3) {
    if (!requireNamespace("mclust", quietly = TRUE)) {
        stop("Package 'mclust' is required for GMM fitting")
    }
    
    # Remove non-finite values
    intensity <- intensity[is.finite(intensity) & intensity > 0]
    
    if (length(intensity) < 100) {
        stop("Insufficient data points for GMM fitting (need at least 100)")
    }
    
    # Fit GMM with error handling
    gmm_fit <- tryCatch({
        mclust::Mclust(intensity, G)
    }, error = function(e) {
        message("Error fitting GMM with G = ", G, ": ", e$message)
        if (G > 1) {
            message("Retrying with G = ", G - 1)
            mclust::Mclust(intensity, G - 1)
        } else {
            stop("Unable to fit GMM even with G = 1")
        }
    })
    
    return(gmm_fit)
}


#' Adjust intensity using GMM-based deconvolution
#'
#' Performs GMM-based background adjustment on intensity data by identifying
#' signal and background components and applying deconvolution using a 
#' correlation-based variance weighting approach.
#'
#' @param intensity Numeric vector of log-transformed intensity values.
#' @param gmm_fit A fitted Mclust object from \code{fitGMMModel}.
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{adj_intensity}: Adjusted intensity values
#'   \item \code{classification}: Component assignments for each pixel
#'   \item \code{posterior}: Posterior probabilities for the signal component
#'   \item \code{G_sig}: Index of the signal component
#'   \item \code{var_weight}: Variance weight (1 - rho) used in adjustment
#'   \item \code{metrics}: Summary metrics for the adjustment
#' }
#'
#' @details
#' The function implements a statistically-principled background correction based on 
#' conditional expectation theory. The method assumes three pixel types:
#' \itemize{
#'   \item Background pixels (C=1): E(U3|X,C=1) = 0
#'   \item Non-specific binding pixels (C=2): E(U3|X,C=2) = 0  
#'   \item Signal pixels (C=3): E(U3|X,C=3) = (mu3 - mu2) + w*(X - mu3)
#' }
#' 
#' The variance weight w = 1 - rho is derived from the correlation structure:
#' \itemize{
#'   \item rho = sign_flip * min(sigma2_bg, sigma2_sig) / sigma2_sig
#'   \item sign_flip = -1 if sigma2_sig < sigma2_bg, else 1
#'   \item This approximates cov(U1+U2, X | C=3) based on bounds and sign constraints
#' }
#' 
#' The final adjustment weights the conditional expectation by the posterior probability
#' P(C=3|X) from the GMM, accounting for uncertainty in component membership.
#' 
#' Note: var_weight can exceed 1 when background variance exceeds signal variance,
#' which amplifies the correction for pixels with high non-specific binding.
#'
#' @export
adjustIntensityGMM <- function(intensity, gmm_fit) {
    G <- gmm_fit$G
    
    # Extract GMM parameters
    comp_means <- gmm_fit$parameters$mean
    comp_proportions <- gmm_fit$parameters$pro
    
    # Handle different variance structures
    if (gmm_fit$parameters$variance$modelName %in% c("E", "EII", "EEI", "EVI")) {
        comp_variance <- rep(gmm_fit$parameters$variance$sigmasq, G)
    } else {
        comp_variance <- gmm_fit$parameters$variance$sigmasq
        if (is.matrix(comp_variance)) {
            comp_variance <- diag(comp_variance)
        }
    }
    
    # Predict component membership
    pred <- predict(gmm_fit, intensity)
    classification <- pred$classification
    posterior_probs <- pred$z
    
    # Determine signal component (highest mean, but check for density overlap)
    xvec <- qnorm(0.99, comp_means[G], sqrt(comp_variance[G]))
    
    if (G > 1 && 
        dnorm(xvec, comp_means[G - 1], sqrt(comp_variance[G - 1])) > 
        dnorm(xvec, comp_means[G], sqrt(comp_variance[G]))) {
        # Second component has higher density at high intensities
        Gbg <- G
        G_sig <- G - 1
    } else {
        Gbg <- G - 1
        G_sig <- G
    }
    
    # Calculate variance weight using correlation coefficient approach
    mu_sig <- comp_means[G_sig]
    mu_bg <- comp_means[Gbg]
    sigma2_sig <- comp_variance[G_sig]
    sigma2_bg <- comp_variance[Gbg]
    
    # Compute covariance approximation
    sign_flip <- ifelse(sigma2_sig < sigma2_bg, -1, 1)
    rho <- sign_flip * min(sigma2_sig, sigma2_bg) / sigma2_sig
    var_weight <- 1 - rho
    
    # Conditional expectation for the signal component
    expected_U_given_signal <- (mu_sig - mu_bg) + var_weight * (intensity - mu_sig)
    
    # Final corrected value: E(U_sig | x) = E(U_sig | x, C=sig) * P(C=sig | x)
    adj_intensity <- expected_U_given_signal * posterior_probs[, G_sig]
    
    # Calculate metrics
    metrics <- list(
        mean_diff = comp_means[G_sig] - comp_means[Gbg],
        var_weight = var_weight,
        neg_adj = sum(adj_intensity < 0, na.rm = TRUE) / length(adj_intensity),
        signal_prop = sum(classification == G_sig) / length(classification)
    )
    
    return(list(
        adj_intensity = adj_intensity,
        classification = classification,
        posterior = posterior_probs[, G_sig],
        G_sig = G_sig,
        var_weight = var_weight,
        metrics = metrics,
        gmm_params = list(
            means = comp_means,
            variances = comp_variance,
            proportions = comp_proportions
        )
    ))
}


#' Apply quantile normalization to adjusted intensities
#'
#' Normalizes adjusted intensities by the 75th percentile of signal pixels.
#'
#' @param adj_intensity Numeric vector of adjusted intensity values.
#' @param classification Integer vector of component classifications.
#' @param G_sig Index of the signal component.
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{normalized}: Quantile-normalized intensities
#'   \item \code{UQ75}: The 75th percentile used for normalization
#' }
#'
#' @importFrom stats quantile
#' @export
quantileNormalizeIntensity <- function(adj_intensity, classification, G_sig) {
    # Get 75th percentile of signal component
    signal_intensities <- adj_intensity[classification == G_sig]
    UQ75 <- quantile(signal_intensities, probs = 0.75, na.rm = TRUE)
    
    # Normalize by UQ75
    normalized <- adj_intensity / UQ75
    
    return(list(
        normalized = normalized,
        UQ75 = UQ75
    ))
}


#' Process rds image channel with GMM normalization
#'
#' Complete pipeline for processing an image channel with GMM-based
#' background normalization.
#'
#' @param file Path to the input RDS file containing intensity matrix.
#' @param channel_name Name of the channel being processed.
#' @param cell_mask Matrix or path to RDS file containing cell segmentation mask.
#' @param cofactor Cofactor for log transformation (default: 150).
#' @param filter_size Median filter size (default: 1).
#' @param bit_depth Bit depth of the input image (default: 16).
#' @param G Number of GMM components (default: 3).
#' @param nsample Number of pixels to sample for GMM fitting (default: 1e5).
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{pixel_data}: Data frame with per-pixel results
#'   \item \code{cell_data}: Data frame with per-cell summaries
#'   \item \code{gmm_summary}: Summary statistics from GMM fitting
#'   \item \code{density_data}: Data for plotting mixture densities
#' }
#'
#' @examples
#' \dontrun{
#' result <- processChannelGMM(
#'   file = "channel_data.rds",
#'   channel_name = "CD45",
#'   cell_mask = "segmentation.rds",
#'   cofactor = 150
#' )
#' }
#'
#' @importFrom EBImage medianFilter
#' @importFrom stats dnorm
#' @export
processChannelGMM <- function(file, channel_name, cell_mask,
                              cofactor = 150, filter_size = 1,
                              bit_depth = 16, G = 3, nsample = 1e5) {
    
    if (!requireNamespace("mclust", quietly = TRUE)) {
        stop("Package 'mclust' is required for GMM processing")
    }
    
    # Load intensity data
    if (is.character(file) && file.exists(file)) {
        intensity_raw <- readRDS(file)
    } else if (is.matrix(file)) {
        intensity_raw <- file
    } else {
        stop("file must be a path to an RDS file or a matrix")
    }
    
    # Load cell mask
    if (is.character(cell_mask) && file.exists(cell_mask)) {
        cell_ids <- readRDS(cell_mask)
        cell_ids <- t(cell_ids)  # Transpose to match orientation
    } else if (is.matrix(cell_mask)) {
        cell_ids <- cell_mask
    } else {
        stop("cell_mask must be a path to an RDS file or a matrix")
    }
    
    # Apply median filter and log transformation
    max_value <- 2^bit_depth - 1
    intensity_filtered <- as.matrix(
        EBImage::medianFilter(intensity_raw / max_value, filter_size)
    ) * max_value
    intensity_log <- log2(intensity_filtered / cofactor + 1)
    
    # Create pixel data frame
    df_pixels <- data.frame(
        channel = channel_name,
        x_coord = rep(seq_len(ncol(intensity_raw)), each = nrow(intensity_raw)),
        y_coord = rep(seq_len(nrow(intensity_raw)), times = ncol(intensity_raw)),
        cell_id = as.vector(cell_ids),
        raw_intensity = as.vector(intensity_raw),
        log_intensity = as.vector(intensity_log)
    )
    
    # Remove non-finite values for GMM fitting
    df_valid <- df_pixels[is.finite(df_pixels$log_intensity) & 
                          df_pixels$log_intensity > 0, ]
    
    # Sample for GMM fitting
    set.seed(123)
    if (nrow(df_valid) > nsample) {
        df_sample <- df_valid[sample(nrow(df_valid), nsample), ]
    } else {
        df_sample <- df_valid
    }
    
    # Fit GMM
    message(paste("Fitting GMM for channel:", channel_name))
    gmm_fit <- fitGMMModel(df_sample$log_intensity, G = G)
    
    # Adjust intensities
    message(paste("Adjusting intensities for channel:", channel_name))
    adjustment <- adjustIntensityGMM(df_pixels$log_intensity, gmm_fit)
    
    df_pixels$adj_intensity <- adjustment$adj_intensity
    df_pixels$classification <- adjustment$classification
    df_pixels$posterior <- adjustment$posterior
    df_pixels$G_sig <- adjustment$G_sig
    
    # Apply quantile normalization
    qn_result <- quantileNormalizeIntensity(
        df_pixels$adj_intensity,
        df_pixels$classification,
        adjustment$G_sig
    )
    
    df_pixels$adj_intensity_qn <- qn_result$normalized
    df_pixels$UQ75 <- qn_result$UQ75
    
    # Aggregate to cell level
    if (!requireNamespace("dplyr", quietly = TRUE)) {
        warning("Package 'dplyr' not available. Skipping cell aggregation.")
        df_cells <- NULL
    } else {
        df_cells <- df_pixels %>%
            dplyr::filter(cell_id != 0) %>%
            dplyr::group_by(cell_id, channel) %>%
            dplyr::summarise(
                x_centroid = mean(x_coord, na.rm = TRUE),
                y_centroid = mean(y_coord, na.rm = TRUE),
                median_raw = median(raw_intensity, na.rm = TRUE),
                median_log = median(log_intensity, na.rm = TRUE),
                median_adj = median(adj_intensity, na.rm = TRUE),
                median_uq = median(adj_intensity_qn, na.rm = TRUE),
                .groups = "drop"
            )
    }
    
    # Create density data for plotting
    G_actual <- gmm_fit$G
    sdev <- sqrt(adjustment$gmm_params$variances)
    if (length(sdev) == 1) {
        sdev <- rep(sdev, G_actual)
    }
    
    intensity_seq <- seq(0, 16, length.out = 500)
    density_components <- sapply(1:G_actual, function(i) {
        dnorm(intensity_seq, 
              adjustment$gmm_params$means[i], 
              sdev[i]) * adjustment$gmm_params$proportions[i]
    })
    
    density_data <- data.frame(
        intensity = intensity_seq,
        mixture_0 = rowSums(density_components),
        density_components
    )
    colnames(density_data)[-1] <- c("0", as.character(1:G_actual))
    
    # GMM summary
    gmm_summary <- data.frame(
        channel = channel_name,
        component = 1:G_actual,
        mean = adjustment$gmm_params$means,
        variance = adjustment$gmm_params$variances,
        proportion = adjustment$gmm_params$proportions,
        G_sig = adjustment$G_sig,
        var_weight = adjustment$var_weight,
        UQ75 = qn_result$UQ75,
        neg_adj = adjustment$metrics$neg_adj,
        signal_prop = adjustment$metrics$signal_prop
    )
    
    message(paste("Completed processing for channel:", channel_name))
    
    return(list(
        pixel_data = df_pixels,
        cell_data = df_cells,
        gmm_summary = gmm_summary,
        density_data = density_data,
        gmm_fit = gmm_fit
    ))
}




#' Create quality control heatmap for GMM results
#'
#' Generates a heatmap visualization of Cohen's d quality control metric
#' across samples and channels.
#'
#' @param gmm_summary_data Data frame containing GMM summary statistics with
#'   columns for sample, channel, and GMM component parameters.
#' @param min_components Minimum number of GMM components required (default: 2).
#' @param samples Optional character vector of sample names to include. If NULL,
#'   uses all samples.
#' @param channels Optional character vector of channel names to include. If NULL,
#'   uses all channels.
#'
#' @return A ggplot2 object representing the QC heatmap.
#'
#' @details
#' This function creates a heatmap showing Cohen's d across samples and channels.
#' Cohen's d provides an interpretable measure of separation between background 
#' and signal components.
#'
#' The function expects a GMM summary data frame with at least the following columns:
#' \itemize{
#'   \item sample: Sample identifier
#'   \item channel: Channel/marker name
#'   \item component: GMM component index
#'   \item mean: Component mean
#'   \item variance: Component variance (or standard deviation)
#' }
#'
#' @examples
#' \dontrun{
#' # Create QC heatmap from GMM results
#' heatmap <- createQCHeatmap(
#'   gmm_summary_data = gmm_summaries
#' )
#' print(heatmap)
#' }
#'
#' @importFrom ggplot2 ggplot aes geom_tile scale_fill_viridis_c theme_minimal
#'   labs theme element_text
#' @importFrom dplyr filter group_by summarise mutate %>%
#' @importFrom tidyr pivot_wider
#' @export
createQCHeatmap <- function(gmm_summary_data, min_components = 2, 
                           samples = NULL, channels = NULL) {
    
    # Check required packages
    check_required_packages(c("ggplot2", "dplyr", "viridis"), "createQCHeatmap")
    
    # Filter samples and channels if specified
    if (!is.null(samples)) {
        gmm_summary_data <- gmm_summary_data %>%
            dplyr::filter(sample %in% samples)
    }
    if (!is.null(channels)) {
        gmm_summary_data <- gmm_summary_data %>%
            dplyr::filter(channel %in% channels)
    }
    
    # Validate GMM data structure
    validate_gmm_summary(gmm_summary_data)
    
    # Calculate metrics for each sample-channel combination
    metrics_df <- gmm_summary_data %>%
        dplyr::group_by(sample, channel) %>%
        dplyr::summarise(
            n_components = dplyr::n(),
            .groups = "drop"
        ) %>%
        dplyr::filter(n_components >= min_components)
    
    if (nrow(metrics_df) == 0) {
        stop("No sample-channel combinations have >= ", min_components, 
             " components", call. = FALSE)
    }
    
    # Join back to get component parameters
    metrics_df <- metrics_df %>%
        dplyr::left_join(gmm_summary_data, by = c("sample", "channel"))
    
    # Calculate Cohen's d between highest two components
    # (typically components 2 vs 3 in a 3-component model,
    # representing non-specific binding vs signal)
    metric_values <- metrics_df %>%
        dplyr::group_by(sample, channel) %>%
        dplyr::summarise(
            max_comp = max(component),
            mu1 = mean[component == max(component) - 1],
            sigma1 = sqrt(variance[component == max(component) - 1]),
            mu2 = mean[component == max(component)],
            sigma2 = sqrt(variance[component == max(component)]),
            .groups = "drop"
        )
    
    # Validate we got the components we need
    if (any(is.na(metric_values$mu1)) || any(is.na(metric_values$mu2))) {
        stop("Could not extract component parameters for all sample-channel pairs. ",
             "Ensure GMM data has at least ", min_components, " components per group.",
             call. = FALSE)
    }
    
    # Compute Cohen's d
    metric_values <- metric_values %>%
        dplyr::mutate(
            metric_value = abs(mu2 - mu1) / sqrt((sigma1^2 + sigma2^2) / 2)
        )
    
    # Create heatmap
    p <- ggplot2::ggplot(metric_values, 
                        ggplot2::aes(x = channel, y = sample, fill = metric_value)) +
        ggplot2::geom_tile(color = "white", linewidth = 0.5) +
        ggplot2::scale_fill_viridis_c(name = "Cohen's d", option = "viridis") +
        ggplot2::theme_minimal() +
        ggplot2::labs(
            title = "QC Heatmap: Cohen's d",
            x = "Channel",
            y = "Sample"
        ) +
        ggplot2::theme(
            axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 10),
            axis.text.y = ggplot2::element_text(size = 10),
            plot.title = ggplot2::element_text(size = 14, face = "bold"),
            legend.position = "right"
        )
    
    return(p)
}
