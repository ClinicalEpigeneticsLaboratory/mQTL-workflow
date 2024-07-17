FROM r-base:4.4.0
ENV DEBIAN_FRONTEND=noninteractive

# Install basic dependencies
RUN apt update && apt upgrade -y
RUN apt install build-essential zlib1g-dev libncurses5-dev libgdbm-dev \
        libnss3-dev libssl-dev libreadline-dev libffi-dev libsqlite3-dev \
        libbz2-dev wget git curl -y

# Download and install htslib-1.9
RUN cd /usr/local/bin && \
    wget https://github.com/samtools/htslib/releases/download/1.9/htslib-1.9.tar.bz2 && \
    tar -vxjf htslib-1.9.tar.bz2 && \
    cd htslib-1.9 && \
    make
ENV PATH="/usr/local/bin/htslib-1.9/:$PATH" 

# Download and install BCFTools-1.9
RUN cd /usr/local/bin && \
    wget https://github.com/samtools/bcftools/releases/download/1.9/bcftools-1.9.tar.bz2 && \
    tar -vxjf bcftools-1.9.tar.bz2 && \
    cd bcftools-1.9 && \
    make
ENV PATH="/usr/local/bin/bcftools-1.9/:$PATH" 

# Download and install Plink 1.9
RUN wget https://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20201019.zip -O /tmp/plink.zip && \
    unzip /tmp/plink.zip -d /usr/local/bin/ && \
    rm /tmp/plink.zip

# Install Illumina CLI 2.1.0
COPY executables/array-analysis-cli-linux-x64-v2.1.0.tar.gz .
RUN gzip -d array-analysis-cli-linux-x64-v2.1.0.tar.gz && \
    tar -xvf array-analysis-cli-linux-x64-v2.1.0.tar && \
    mv array-analysis-cli-linux-x64-v2.1.0/ /usr/local/bin/ && \
    rm array-analysis-cli-linux-x64-v2.1.0.tar 
ENV PATH="/usr/local/bin/array-analysis-cli-linux-x64-v2.1.0/array-analysis-cli:$PATH" 

# Install VEP
RUN apt install perl libdbi-perl libdbd-mysql-perl liblist-moreutils-perl libwww-perl libmodule-build-perl -y
RUN cd /usr/local/bin && git clone https://github.com/Ensembl/ensembl-vep.git && \
    cd ensembl-vep && \
    perl INSTALL.pl --AUTO a
ENV PATH="/usr/local/bin/ensembl-vep/:$PATH" 


# # Set the working directory
WORKDIR /workspace

# # Copy the codebase into the container.
COPY requirements.R poetry.lock pyproject.toml info.sh .

# Install R components
RUN apt install -y libcurl4-openssl-dev libssl-dev libxml2-dev
RUN Rscript requirements.R

# Install Python
RUN wget https://www.python.org/ftp/python/3.10.0/Python-3.10.0.tgz && \
    tar -xvf Python-3.10.0.tgz && \
    cd Python-3.10.0 && \
    ./configure --enable-optimizations && \ 
    make -j 4 && \
    make altinstall

# # Install requirements
RUN python3.10 -m pip install poetry && python3.10 -m poetry export --without-hashes --format=requirements.txt > requirements.txt
RUN python3.10 -m pip install -r requirements.txt

# Run a command by default
CMD ["sh", "info.sh"]
