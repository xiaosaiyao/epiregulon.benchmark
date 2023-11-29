#' @import basilisk
venv1 <- BasiliskEnvironment(envname="cellOracle",
                              pkgname="epiregulon.benchmark",
                              packages=c("bedtools=2.30.0", "numpy==1.23.5", "cython==0.29.33"),
                                         channels = c("conda-forge","pypi",
                                                      "bioconda"),
                            pip=c("velocyto==0.17.17", "celloracle==0.12.1"))

setBasiliskCheckVersions(FALSE) # allow for pip installation from github
venv2 <- BasiliskEnvironment(envname="scenicplus",
                             pkgname="epiregulon.benchmark",
                             pip=c("git+https://github.com/aertslab/scenicplus", "scikit-learn==1.3.1",
                                   "pyranges==0.0.127","matplotlib==3.7.3",
                                   "seaborn==0.12.2", "attr==0.3.2", "ray==2.7.0","tqdm==4.66.1"
                                   ),
                             packages = c("python=3.8", "scanpy=1.9.5", "cython=0.29.36", "numpy=1.24.4",
                                          "pandas=1.5.0", "scipy=1.10.1", "pybiomart=0.2.0"),
                             channels = c("conda-forge", "bioconda", "anaconda"))
