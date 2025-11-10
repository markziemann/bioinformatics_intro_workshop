FROM bioconductor/bioconductor_docker:3.21-R-4.5.2

# Update apt-get
RUN apt-get update \
        && apt-get upgrade -y \
        && apt-get install -y nano git  libncurses-dev \
        ## Remove packages in '/var/cache/' and 'var/lib' to remove side-effects of apt-get update
        && apt-get clean \
        && rm -rf /var/lib/apt/lists/*

# Install CRAN packages
RUN Rscript -e 'install.packages(c("eulerr","gplots","kableExtra","parallel","RColorBrewer","stringi","HGNChelper","beeswarm","gridExtra","png","RhpcBLASctl","tictoc"))'

# Install bioconductor packages
RUN Rscript -e 'BiocManager::install(c("biomaRt","clusterProfiler","DESeq2","edgeR","fgsea","getDEE2","limma","mitch"))'

RUN Rscript -e 'install.packages(c("vioplot","Hmisc","rmdformats"))'

RUN Rscript -e 'install.packages("caret", dependencies = c("Depends", "Suggests"))'

# get a clone of the codes
#RUN git clone https://github.com/markziemann/myproject.git
RUN mkdir /myproject

# Set the container working directory
ENV DIRPATH /myproject
WORKDIR $DIRPATH
