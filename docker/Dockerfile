FROM rocker/drd

RUN apt-get update && apt-get install \
	-y --no-install-recommends \
	libssl-dev
RUN apt-get update && apt-get -y --no-install-recommends install libcurl3-gnutls git
 
RUN Rscript -e 'utils::install.packages(c("rmarkdown", "knitr", "coda", "mvtnorm", "tinytest", "mcmc", "covr"))'

RUN apt-get update && \
	apt install -y --no-install-recommends --allow-change-held-packages libicu67=67.1-7 && \
	apt-get install -y --no-install-recommends -f libxml2-dev && \
	Rscript -e 'utils::install.packages(c("xml2"))'


CMD ["bash"]

