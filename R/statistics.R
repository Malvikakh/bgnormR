#' Calculate Cohen's d effect size
#'
#' Computes Cohen's d, a standardized measure of effect size between two groups.
#'
#' @param mu1 Mean of the first group.
#' @param mu2 Mean of the second group.
#' @param sigma1 Standard deviation of the first group.
#' @param sigma2 Standard deviation of the second group.
#'
#' @return Numeric value representing Cohen's d effect size.
#'
#' @details
#' Cohen's d is calculated as the difference between two means divided by the
#' pooled standard deviation:
#' \deqn{d = \frac{|\mu_2 - \mu_1|}{\sqrt{(\sigma_1^2 + \sigma_2^2) / 2}}}
#'
#' Interpretation guidelines:
#' \itemize{
#'   \item Small effect: |d| = 0.2
#'   \item Medium effect: |d| = 0.5
#'   \item Large effect: |d| = 0.8
#' }
#'
#' @examples
#' # Calculate Cohen's d between two normal distributions
#' cohens_d(mu1 = 5, mu2 = 7, sigma1 = 1.5, sigma2 = 2.0)
#'
#' @export
cohens_d <- function(mu1, mu2, sigma1, sigma2) {
    inputs <- list(mu1 = mu1, mu2 = mu2, sigma1 = sigma1, sigma2 = sigma2)
    if (!all(sapply(inputs, is.numeric))) {
        stop("All inputs must be numeric")
    }
    if (sigma1 <= 0 || sigma2 <= 0) {
        stop("Standard deviations must be positive")
    }
    
    pooled_sd <- sqrt((sigma1^2 + sigma2^2) / 2)
    d <- abs(mu2 - mu1) / pooled_sd
    
    return(d)
}


#' Calculate weighted Cohen's d
#'
#' Computes a variance-weighted version of Cohen's d effect size.
#'
#' @param mu1 Mean of the first group.
#' @param mu2 Mean of the second group.
#' @param sigma1 Standard deviation of the first group.
#' @param sigma2 Standard deviation of the second group.
#' @param var_weight Variance weight to apply.
#'
#' @return Numeric value representing weighted Cohen's d.
#'
#' @details
#' This is a weighted variant of Cohen's d that accounts for differences in
#' variance structure between groups:
#' \deqn{d_w = \frac{|\mu_2 - \mu_1| \times w}{\sqrt{(\sigma_1^2 + \sigma_2^2) / 2}}}
#' where w is the variance weight.
#'
#' @examples
#' weighted_cohens_d(mu1 = 5, mu2 = 7, sigma1 = 1.5, sigma2 = 2.0, var_weight = 1.2)
#'
#' @export
weighted_cohens_d <- function(mu1, mu2, sigma1, sigma2, var_weight) {
    inputs <- list(mu1 = mu1, mu2 = mu2, sigma1 = sigma1, sigma2 = sigma2, var_weight = var_weight)
    if (!all(sapply(inputs, is.numeric))) {
        stop("All inputs must be numeric")
    }
    if (sigma1 <= 0 || sigma2 <= 0) {
        stop("Standard deviations must be positive")
    }
    
    pooled_sd <- sqrt((sigma1^2 + sigma2^2) / 2)
    d_weighted <- abs(mu2 - mu1) * var_weight / pooled_sd
    
    return(d_weighted)
}
