#' @importFrom basilisk BasiliskEnvironment
#' @importFrom basilisk basiliskStart basiliskRun basiliskStop
env1 <- BasiliskEnvironment("env1", pkgname="APL",
    packages=c("numpy==1.19.5", "pytorch==1.10.1"))
