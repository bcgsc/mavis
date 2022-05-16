FROM python:3.7-slim-buster

WORKDIR /app

RUN apt-get update && \
    apt-get upgrade -y && \
    apt-get install -y git wget make gcc libz-dev

# pysam dependencies
RUN apt-get install -y libncurses5-dev zlib1g-dev libbz2-dev libncursesw5-dev liblzma-dev

# install BWA
RUN git clone https://github.com/lh3/bwa.git && \
    cd bwa && \
    git checkout v0.7.17 && \
    make && \
    cd .. && \
    mv bwa/bwa /usr/local/bin

# install minimap2
RUN git clone https://github.com/lh3/minimap2.git && \
    cd minimap2 && \
    git checkout v2.24 && \
    make && \
    cd .. && \
    mv minimap2/minimap2.1 /usr/local/bin

# install blat dependencies
RUN apt-get install -y libcurl4

# install blat
RUN wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/blat/blat && \
    chmod a+x blat && \
    mv blat /usr/local/bin

# install wtdbg2
RUN git clone https://github.com/ruanjue/wtdbg2.git && \
    cd wtdbg2 && \
    make && \
    cd .. && \
    mv wtdbg2/wtdbg2 /usr/local/bin

COPY setup.py setup.py
COPY setup.cfg setup.cfg
COPY MANIFEST.in MANIFEST.in
COPY pyproject.toml pyproject.toml
COPY src src
COPY LICENSE LICENSE
COPY README.md README.md

# install python package
RUN pip install -U setuptools pip wheel
RUN pip install .
RUN which mavis
ENTRYPOINT [ "mavis" ]
