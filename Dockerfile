FROM rocker/verse:4.2.3

# define input and output directories for mounting
ENV INPUT_DIR=/input
ENV OUTPUT_DIR=/output

RUN R -e 'install.packages("remotes")'
RUN R -e 'remotes::install_cran("rlang")'
RUN apt-get update && apt-get install -y libcurl4-openssl-dev libhdf5-dev libicu-dev libssl-dev make pandoc zlib1g-dev && rm -rf /var/lib/apt/lists/*
RUN mkdir -p /usr/local/lib/R/etc/ /usr/lib/R/etc/
RUN echo "options(repos = c(CRAN = 'https://cran.rstudio.com/'), download.file.method = 'libcurl', Ncpus = 4)" | tee /usr/local/lib/R/etc/Rprofile.site | tee /usr/lib/R/etc/Rprofile.site
RUN Rscript -e 'remotes::install_version("rbokeh")'

RUN apt-get update && apt-get install patch
RUN Rscript -e "remotes::install_github('adnaniazi/tailfindr', ref='master')"

# Install VBZ plugin
RUN apt-get update && apt-get install -y wget
RUN wget https://github.com/nanoporetech/vbz_compression/releases/download/v1.0.1/ont-vbz-hdf-plugin-1.0.1-Linux-x86_64.tar.gz && \
    tar -xzf ont-vbz-hdf-plugin-1.0.1-Linux-x86_64.tar.gz -C /usr/local/
ENV HDF5_PLUGIN_PATH "/usr/local/ont-vbz-hdf-plugin-1.0.1-Linux/usr/local/hdf5/lib/plugin"

CMD R -e 'library(tailfindr)'
