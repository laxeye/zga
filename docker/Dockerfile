# Define base image
FROM continuumio/miniconda3

MAINTAINER laxeye

# Create Conda environment from the YAML file
COPY environment.yml .
RUN conda env create -f environment.yml
RUN conda clean -ay

# Download light bakta database
# https://github.com/oschwengers/bakta?tab=readme-ov-file#database-download
RUN bakta_db download --output ${HOME}/bakta --type light

# Override default shell and use bash
SHELL ["conda", "run", "-n", "zga", "/bin/bash", "-c"]

# Python program to run in the container
ENTRYPOINT ["conda", "run", "-n", "zga", "zga"]
