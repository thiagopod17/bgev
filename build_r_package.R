# building the r package
rm(list=ls())
code_files = c("package/bgev_functions.R","package/dist_check.R")
package_name = 'bgev'
path_for_package_files = 'bgev_cran'
package.skeleton(name=package_name, 
                 path=path_for_package_files, 
                 code_files = code_files)
