#' @import basilisk
venv1 <- BasiliskEnvironment(envname="cellOracle",
                              pkgname="epiregulon.benchmark",
                              packages=c("bedtools=2.30.0", "numpy==1.23.5", "cython==0.29.33"),
                                         channels = c("conda-forge","pypi",
                                                      "bioconda"),
                            pip=c("velocyto==0.17.17", "celloracle==0.12.1")
)


venv2 <- BasiliskEnvironment(envname="scenicplus",
                            pkgname="epiregulon.benchmark",
                            pip=c("ray==2.5.0"),
                            packages = "python=3.8",
                            channels = "conda-forge",
                            paths = "scenicplus"
)
