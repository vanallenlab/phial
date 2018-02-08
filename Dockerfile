FROM vanallenlab/r-base:3.4.2

RUN Rscript -e 'install.packages("optparse", repos="http://ftp.ussg.iu.edu/CRAN/")'
RUN Rscript -e 'install.packages("gplots", repos="http://ftp.ussg.iu.edu/CRAN/")'
RUN Rscript -e 'install.packages("Nozzle.R1", repos="http://ftp.ussg.iu.edu/CRAN/")'

WORKDIR /

COPY . /
