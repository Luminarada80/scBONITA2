# Use the official Miniconda3 base image
FROM continuumio/miniconda3

# Set the working directory inside the container
WORKDIR /app

# Copy the environment.yml file to the working directory
COPY environment.yml /app/environment.yml

# Install mamba and create the Conda environment
RUN conda install mamba -n base -c conda-forge \
    && mamba env create -f /app/environment.yml \
    && conda clean -afy

# Initialize conda in the shell
RUN echo "source /opt/conda/etc/profile.d/conda.sh" >> ~/.bashrc
RUN echo "conda activate scBonita" >> ~/.bashrc

# Copy the scBONITA code into the container
COPY . /app

# Set the working directory inside the container
WORKDIR /app

# Run the default command inside the activated conda environment
CMD ["bash", "-c", "source /opt/conda/etc/profile.d/conda.sh && conda activate scBonita && bash bash_scripts/local_george_hiv.sh"]
