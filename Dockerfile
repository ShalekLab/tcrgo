FROM debian:latest
# ORIGINAL MAINTAINER Conda Development Team <conda@continuum.io>

RUN apt-get -qq update && apt-get -qq -y install curl unzip bzip2 ca-certificates default-jre libgomp1 build-essential zlib1g-dev \
	&& curl -L -s -o "/Drop-seq_tools-2.4.0.zip" "https://github.com/broadinstitute/Drop-seq/releases/download/v2.4.0/Drop-seq_tools-2.4.0.zip" \
	&& unzip /Drop-seq_tools-2.4.0.zip \
	&& rm /Drop-seq_tools-2.4.0.zip \
	&& mkdir /software \
	&& mv /Drop-seq_tools-2.4.0 /software/dropseq \
	&& chmod a+rx /software/dropseq/* \
	&& curl -sSL https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -o /tmp/miniconda.sh \
	&& bash /tmp/miniconda.sh -bfp /usr/local \
	&& rm -rf /tmp/miniconda.sh \
	&& conda install -y python=3.8.5 \
	&& conda update conda \
	&& apt-get -qq -y remove curl bzip2 \
	&& apt-get -qq -y autoremove \
	&& apt-get autoclean \
	&& rm -rf /var/lib/apt/lists/* /var/log/dpkg.log \
	&& conda clean --all --yes \
	&& conda config --add channels conda-forge --add channels bioconda \
	&& conda install -y samtools=1.3.1 bowtie2=2.4.1 pysam=0.16.0.1 pandas=1.1.2

COPY transform_read_data.sh /scripts/
COPY *.py /scripts/
COPY tcrgo /scripts/tcrgo