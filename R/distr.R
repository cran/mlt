
.Normal <- function()
    list(p = pnorm, d = dnorm, q = qnorm, 
         ### see also MiscTools::ddnorm
         dd = function(x) -dnorm(x = x) * x,
         ddd = function(x) dnorm(x = x) * (x^2 - 1), 
         name = "normal")

.Logistic <- function()
    list(p = plogis, d = dlogis, q = qlogis,
         dd = function(x) {
             ex <- exp(x)
             (ex - exp(2 * x)) / (1 + ex)^3
         },
         ddd = function(x) {
             ex <- exp(x)
             (ex - 4*(exp(2 * x)) + exp(3 * x)) / (1 + ex)^4
         },
         name = "logistic")

.MinExtrVal <- function()
    list(p = function(x) 1 - exp(-exp(x)),
         q = function(p) log(-log(1 - p)),
         d = function(x, log = FALSE) {
             ret <- x - exp(x)
             if (!log) return(exp(ret))
             ret
         },
         dd = function(x) {
             ex <- exp(x)
             (ex - exp(2 * x)) / exp(ex)
         },
         ddd = function(x) {
             ex <- exp(x)
             (ex - 3*exp(2 * x) + exp(3 * x)) / exp(ex)
         },
         name = "minimum extreme value")

.distr <- function(which = c("Normal", "Logistic", 
                             "MinExtrVal")) {
    which <- match.arg(which)
    do.call(paste(".", which, sep = ""), list())
}
