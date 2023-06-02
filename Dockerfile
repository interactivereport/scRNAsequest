FROM debian:stable-20220316

USER root
ENV TZ=America/New_York
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone

# Install miniconda
ENV PATH=$PATH:/root/miniconda3/bin
ARG PATH=$PATH:/root/miniconda3/bin

RUN apt-get update

RUN apt-get install -y wget && rm -rf /var/lib/apt/lists/*

RUN wget \
    https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh --no-check-certificate \
    && mkdir /root/.conda \
    && bash Miniconda3-latest-Linux-x86_64.sh -b \
    && rm -f Miniconda3-latest-Linux-x86_64.sh 

RUN conda --version

# Copy source files of scRNASequest
RUN mkdir /home/scRNASequest
WORKDIR /home/scRNASequest
COPY . .

# Install scRNASequest env
RUN bash install.sh

ENV PATH=$PATH:/home/scRNASequest

# Keep the container running
CMD while true; do sleep 1000; done
