url<-"https://cran.r-project.org/src/contrib/Archive/SNPfiltR/SNPfiltR_1.0.1.tar.gz"
pkgFile<- "SNPfiltR_1.0.1.tar.gz"
download.file(url = url, destfile = pkgFile)

#I unzipped file in my Windows console

#install package
install.packages(c("Rtsne", "ggridges")) #install dependencies
install.packages(pkgs = pkgFile, type = "source", repos=NULL)
