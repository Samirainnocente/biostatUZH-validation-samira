
#' Quantile of a Kaplan-Meier estimate of a survival function
#' 
#' In survival analysis, often quantiles (e.g., the median) of an
#' estimated survival function are of interest. The \code{\link[survival]{survfit}}
#' function in the  'survival' packages of version \eqn{< 2.35} computes median time to event.
#' However, the output of that function is such that the median is not
#' conveniently accessible by the user. This function
#' makes the median and any other quantile accessible by the user.
#' 
#' @param time Event times, censored or observed.
#' @param event Censoring indicator, 1 for event, 0 for censored.
#' @param group Indicates groups to compute individual quantile for each group.
#' Default is \code{NULL}.
#' @param quantile Quantile to be computed. Real number in \eqn{[0, 1]}.
#' @param conf.level Significance level for confidence interval for the time to
#' event quantile.
#' @param conf.type Type of confidence interval to be computed. For possible
#' choices see above, and for specifications regarding the different options
#' see \code{\link[survival]{survfit}}.
#' @return
#' A list containing a matrix with columns:
#' \item{n}{Number of observations used.}
#' \item{events}{Number of events.}
#' \item{quantile}{Quantile estimate.}
#' \item{lower.ci}{Lower limit of confidence interval.}
#' \item{upper.ci}{Upper limit of confidence interval.}\cr
#' and the p-value of the test for the difference between the survival curves;
#' see \code{\link[survival]{survdiff}} for more information.
#' @author Kaspar Rufibach
#' @note The function is based on the function \code{\link[survival]{survfit}}.
#' @seealso \code{\link[survival]{survfit}}, \code{\link[survival]{survdiff}}
#' @references Confidence intervals are computed according to
#' 
#' Brookmeyer, R. and Crowley, J. (1982).  A Confidence Interval for the Median
#' Survival Time.  \emph{Biometrics}, \bold{38}, 29--41.
#' @keywords htest survival
#' @examples
#' 
#' ## use Acute Myelogenous Leukemia survival data contained in package 'survival'
#' time <- leukemia[, 1]
#' event <- leukemia[, 2]
#' group <- as.factor(leukemia[, 3])
#' plot(survfit(Surv(time, event) ~ group, conf.type = "none"), mark = "/", col = 1:2)
#' 
#' ## median time to event
#' quantileKM(time, event, group = group, quantile = 0.5, conf.level = 0.95, 
#'            conf.type = "log")
#'  
#' ## comparison to standard function (median time to event not accessible by user)
#' survfit(Surv(time, event) ~ group, conf.int = 0.95)
#' 
#' ## compute 0.25 quantile
#' quantileKM(time = time, event = event, group = group, quantile = 0.25,
#'            conf.type = "log")
#'
#' @import survival
#' @export
quantileKM <- function(time, event, group = NULL, quantile = 0.5, conf.level = 0.95,
                       conf.type = c("log-log", "log", "plain","none")){
  
  stopifnot(is.numeric(time),
            length(time) > 0,
            is.numeric(event),
            length(event) == length(time))
  if(!is.null(group)){
    group <- as.factor(group)
    stopifnot(is.finite(group),
              length(group) == length(time))
  }
  stopifnot(is.numeric(quantile), length(quantile) == 1,
            is.finite(quantile),
            0 <= quantile, quantile <= 1,
            is.numeric(conf.level), length(conf.level) == 1,
            is.finite(conf.level),
            0 < conf.level, conf.level < 1)
  
  stopifnot(!is.null(conf.type))
  conf.type <- match.arg(conf.type)
  
  s.obj <- Surv(time, event)
  
  if (!is.null(group)) {
    group.f <- as.factor(group)
    n.level <- length(levels(group.f))
    quant.mat <- matrix(NA, ncol = 5, nrow = n.level)
    sdiff <- survdiff(s.obj ~ group.f)
    
    ## p-value log-rank test
    p.val <- 1 - pchisq(sdiff$chisq, df = n.level - 1)
    quant.surv <- survfit(s.obj ~ group.f, conf.int = conf.level,
                          conf.type = conf.type)
    
    ## quantile of event time, incl. confidence intervals
    for (j in 1:n.level){
      tmp <- summary(quant.surv[j])
      quant.mat[j, 1] <- quant.surv[j]$n
      quant.mat[j, 2] <- sum(quant.surv[j]$n.event)
      
      ## add Inf to omit warnings in case quantile is not reached by 
      ## survival curve (or pointwise ci curve)
      quant.mat[j, 3] <- min(Inf, tmp$time[tmp$surv <= quantile], na.rm = TRUE)
      quant.mat[j, 4] <- min(Inf, tmp$time[tmp$lower <= quantile], na.rm = TRUE)
      quant.mat[j, 5] <- min(Inf, tmp$time[tmp$upper <= quantile], na.rm = TRUE)
    }        
    
    dimnames(quant.mat)[[1]] <- paste("group = ", levels(group.f), sep = "")
  }
  
  if (is.null(group)) {
    p.val <- NA
    n.level <- 1
    quant.mat <- matrix(NA, ncol = 5, nrow = n.level)
    quant.surv <- survfit(s.obj ~ 1, conf.int = conf.level, conf.type = conf.type)
    tmp <- quant.surv
    quant.mat[1, 1] <- tmp$n
    quant.mat[1, 2] <- sum(tmp$n.event)
    
    ## add Inf to omit warnings in case quantile is not reached by 
    ## survival curve
    quant.mat[1, 3] <- min(Inf, tmp$time[tmp$surv <= quantile], na.rm = TRUE)
    quant.mat[1, 4] <- min(Inf, tmp$time[tmp$lower <= quantile], na.rm = TRUE)
    quant.mat[1, 5] <- min(Inf, tmp$time[tmp$upper <= quantile], na.rm = TRUE)        
  }    
  
  dimnames(quant.mat)[[2]] <- c("n", "events", "quantile", "lower.ci", "upper.ci")      
  list("quantities" = quant.mat, "p.val" = p.val)
}
