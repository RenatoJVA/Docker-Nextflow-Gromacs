FROM nvidia/cuda:12.9.1-cudnn-devel-ubuntu24.04

# =============================
# Variables de entorno
# =============================
ENV DEBIAN_FRONTEND=noninteractive \
  PATH=/root/anaconda3/bin:$PATH \
  GMX_MAXBACKUP=999

# =============================
# Instalar dependencias
# =============================
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
  ca-certificates \
  && rm -rf /var/lib/apt/lists/*

# =============================
# Instalar Nextflow
# =============================
WORKDIR /usr/local/bin
RUN curl -s https://get.nextflow.io | bash && \
  chmod +x nextflow

WORKDIR /opt

# =============================
# Descargar GROMACS
# =============================
RUN wget https://ftp.gromacs.org/gromacs/gromacs-2024.4.tar.gz && \
  tar xfz gromacs-2024.4.tar.gz && \
  rm gromacs-2024.4.tar.gz

# =============================
# Compilar versión SINGLE (SP)
# =============================
RUN cd gromacs-2024.4 && mkdir build && cd build && \
  CC=gcc CXX=g++ cmake .. \
  -DGMX_OPENMP=ON \
  -DGMX_BUILD_OWN_FFTW=ON \
  -DREGRESSIONTEST_DOWNLOAD=ON \
  -DGMX_GPU=CUDA \
  -DCUDA_TOOLKIT_ROOT_DIR=/usr/local/cuda \
  -DGMX_FFT_LIBRARY=fftw3 \
  -DGMX_PYTHON_PACKAGE=OFF \
  -DGMX_PREFER_STATIC_LIBS=ON && \
  make -j"$(nproc)" && make install

# =============================
# Compilar versión DOUBLE (DP)
# =============================
RUN cd gromacs-2024.4 && rm -rf build && mkdir build && cd build && \
  CC=gcc CXX=g++ cmake .. \
  -DGMX_OPENMP=ON \
  -DGMX_BUILD_OWN_FFTW=ON \
  -DREGRESSIONTEST_DOWNLOAD=ON \
  -DGMX_FFT_LIBRARY=fftw3 \
  -DGMX_PYTHON_PACKAGE=OFF \
  -DGMX_PREFER_STATIC_LIBS=ON \
  -DGMX_DOUBLE=on && \
  make -j"$(nproc)" && make install

# =============================
# Configurar entorno
# =============================
RUN echo -e "\n#########################################################################\n########################### Gromacs #####################################\nsource /usr/local/gromacs/bin/GMXRC\nexport GMX_MAXBACKUP=999" \
  >> /etc/bash.bashrc

ENV PATH=/usr/local/gromacs-sp/bin:/usr/local/gromacs-dp/bin:$PATH

# =============================
# Carpeta de trabajo
# =============================
WORKDIR /workspace

# =============================
# Entrada por defecto
# =============================
CMD ["bash"]% 