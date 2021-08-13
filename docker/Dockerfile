FROM rocker/drd

RUN apt-get update && apt-get install \
	-y --no-install-recommends \
	libssl-dev
RUN apt-get update && apt-get -y --no-install-recommends install libcurl3-gnutls git
RUN Rscript -e 'utils::install.packages(c("rmarkdown", "knitr", "coda", "mvtnorm", "tinytest", "mcmc", "covr"))'

CMD ["bash"]

