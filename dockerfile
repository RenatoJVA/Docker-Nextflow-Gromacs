FROM nvidia/cuda:12.4.1-base-ubuntu22.04

# Variables
ENV DEBIAN_FRONTEND=noninteractive \
    PATH=/root/anaconda3/bin:$PATH \
    GMX_MAXBACKUP=999

# Instalar dependencias
RUN apt-get update && apt-get install -y \
    build-essential \
    cmake \
    git \
    wget \
    curl \
    vim \
    tree \
    btop \
    bzip2 \
    tmux \
    openjdk-17-jre \
    libfftw3-dev \
    mpich \
    && rm -rf /var/lib/apt/lists/*

# Instalar Nextflow
WORKDIR /usr/local/bin
RUN curl -s https://get.nextflow.io | bash && \
    chmod +x nextflow

# Instalar Anaconda
RUN wget https://repo.anaconda.com/archive/Anaconda3-2024.10-1-Linux-x86_64.sh -O /tmp/anaconda.sh && \
    bash /tmp/anaconda.sh -b -p /root/anaconda3 && \
    rm /tmp/anaconda.sh && \
    /root/anaconda3/bin/conda init --all

# Instalar GROMACS 2024.4 (doble precisión)
WORKDIR /opt
RUN wget https://ftp.gromacs.org/gromacs/gromacs-2024.4.tar.gz && \
    tar xfz gromacs-2024.4.tar.gz && \
    cd gromacs-2024.4 && \
    mkdir build && cd build && \
    CC=gcc CXX=g++ cmake .. \
      -DGMX_OPENMP=ON \
      -DGMX_BUILD_OWN_FFTW=ON \
      -DREGRESSIONTEST_DOWNLOAD=ON \
      -DGMX_GPU=CUDA \
      -DCUDA_TOOLKIT_ROOT_DIR=/usr/local/cuda \
      -DGMX_FFT_LIBRARY=fftw3 \
      -DGMX_PYTHON_PACKAGE=OFF \
      -DGMX_PREFER_STATIC_LIBS=ON \
      -DGMX_DOUBLE=ON && \
    make -j$(nproc) && \
    make -j$(nproc) check && \
    make install && \
    cd /opt && rm -rf gromacs-2024.4 gromacs-2024.4.tar.gz

# Añadir GROMACS a entorno global
RUN echo "source /usr/local/gromacs/bin/GMXRC" >> /etc/bash.bashrc

# Usar bash como shell para que cargue conda
SHELL ["/bin/bash", "-c"]

# Carpeta de trabajo por defecto
WORKDIR /workspace

# Verificar instalación
RUN java -version && nextflow -version && gmx --version && conda --version

# Entrada por defecto
CMD ["/bin/bash"]
