args <- commandArgs(trailingOnly = T)
cat(sprintf("Args: %s\n", args))
setwd(args[4])
Sweave(args[1])
