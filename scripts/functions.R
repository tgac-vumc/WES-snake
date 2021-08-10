
#Author: Matias Mendeville
#######################
### Custom R functions
#######################


# To source all (.R-) files in a directory
sourceDir <- function(path, trace = TRUE, ...) {
         for (nm in list.files(path, pattern = "[.][RrSsQq]$")) {
            if(trace) cat(nm,":")
            source(file.path(path, nm), ...)
            if(trace) cat("\n")
         }
      }
