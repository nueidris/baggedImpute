#' @import stats
#' @importFrom imputeTS na_interpolation na_kalman na_ma na_seadec
#' @export

baggedImpute <- function(data, bootrep = 1000, blocklen = 2, BBM = c("MBB", "NBB", "CBB"), misstr = c("ams", "ivs"), interpolation) {
    bbm <- data.frame(
        "bbm" = c("MBB", "NBB", "CBB"),
        "formulas" = c("c(GenCBB)", "c(GenMBB)", "c(GenNBB)")
    )
    methods <- data.frame(
      "methods" = c(
        "auto.arima", "StructTS", "linear_i",
        "spline_i", "stine_i", "simple_ma", "linear_ma",
        "exponential_ma", "seadec_kalman", "seadec_ma",
        "seadec_random"
      ),
      "formulas" = c(
        "na_kalman(boot.result, model = 'auto.arima', smooth = TRUE, nit = -1)",
        "na_kalman(boot.result, model = 'StructTS', smooth = TRUE, nit = -1)",
        "na_interpolation(boot.result, option = 'linear')",
        "na_interpolation(boot.result, option = 'spline')",
        "na_interpolation(boot.result, option = 'stine')",
        "na_ma(boot.result, k=3, weighting = 'simple')",
        "na_ma(boot.result, k=3, weighting = 'linear')",
        "na_ma(boot.result, k=3, weighting = 'exponential')",
        "na_seadec(boot.result, algorithm = 'kalman')",
        "na_seadec(boot.result, algorithm = 'ma')",
        "na_seadec(boot.result, algorithm = 'random')"
      )
    )

    if (BBM %in% bbm$bbm) {
        BootFunc <- eval(parse(text = bbm$formulas[bbm$bbm == BBM]))
    } else stop("Error: BBM not available!")
    
    boot.result <- BootFunc[[1]](Xt = data, l = blocklen, r = bootrep, irreg = misstr)
    if (interpolation %in% methods$methods) {
      z <- eval(parse(text = methods$formulas[methods$methods == interpolation]))
    } else stop("Error: method not available!")

    listNA <- which(is.na(data))
    imputedSeries <- data
    imputedSeries[listNA] <- (apply(z , 1, mean , na.rm = TRUE))[listNA]
    return(imputedSeries)
}
