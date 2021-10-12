# Define base image
FROM continuumio/miniconda3

MAINTAINER laxeye

# Create Conda environment from the YAML file
COPY environment.yml .
RUN conda env create -f environment.yml
RUN conda clean -ay

# Download more DFAST databases and proteins
# You should patch default_config.py to use them
# RUN dfast_file_downloader.py --cdd Cdd_NCBI Prk Smart--hmm Pfam dbCAN
# RUN dfast_file_downloader.py --protein bifido cyanobase ecoli lab

# Override default shell and use bash
SHELL ["conda", "run", "-n", "zga", "/bin/bash", "-c"]

# Python program to run in the container
ENTRYPOINT ["conda", "run", "-n", "zga", "zga"]
