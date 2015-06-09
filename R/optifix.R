#' Fixed parameter optimisation
#'
#' Optifix is a function originally written by Barry Rowlingson (b.rowlingsno at lancaster.ac.uk) that acts as a wrapper for the optim function, allowing the user to fix any specified parameters during optimisation.
#' @param par a vector of the parameters to be optimised
#' @param fixed a vector of boolean values corresponding to the values of par. If the corresponding index is TRUE, then the parameter in par will be fixed during optimisation.
#' @param fn the cost function to be optimised over, as in optim
#' @param gr as in optim
#' @return returns as with optim, but with a few extra variables to summarise the fixed parameters
#' @seealso \code{\link{optim}}
#' @export
#' @examples
#' fit <- optifix(c(1,1),c(TRUE,FALSE),cost_function,data=example_data)
optifix <- function(par, fixed, fn, gr = NULL, ...,
           method = c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN"),
           lower = -Inf, upper = Inf,
           control = list(), hessian = FALSE){
   force(fn)
   force(fixed)
   .npar=length(par)
   .fixValues = par[fixed]

   .parStart = par[!fixed]
   
   .fn <- function(par,...){
     .par = rep(NA,sum(!fixed))
     .par[!fixed] = par
     .par[fixed] = .fixValues
     fn(.par,...)
   }

   if(!is.null(gr)){
     .gr <- function(par,...){
       .gpar = rep(NA,sum(!fixed))
       .gpar[!fixed] = par
       .gpar[fixed] = .fixValues
       gr(.gpar,...)[!fixed]
     }
   }else{
     .gr <- NULL
   }

   .opt = optim(.parStart,.fn,.gr,...,method=method,lower=lower,control=control,hessian=hessian)

   .opt$fullpars = rep(NA,sum(!fixed))
   .opt$fullpars[fixed]=.fixValues
   .opt$fullpars[!fixed]=.opt$par
   .opt$fixed = fixed
   return(.opt)
   
 }
