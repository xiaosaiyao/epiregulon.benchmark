venv <- BasiliskEnvironment(envname="cellOracle",
                              pkgname="epiregulon.benchmark",
                              packages=c("bedtools=2.30.0", "numpy==1.23.5", "cython==0.29.33"),
                                         channels = c("conda-forge","pypi",
                                                      "bioconda"),
                            pip=c("velocyto==0.17.17", "celloracle==0.12.1")
)


