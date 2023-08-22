FROM r-base:latest

# Install system dependencies
RUN apt-get update && \
    apt-get install -y libcurl4-openssl-dev libcurl4-doc libidn-dev libkrb5-dev libldap2-dev librtmp-dev libssh2-1-dev libssl-dev libpoppler-cpp-dev && \
    rm -rf /var/lib/apt/lists/*

# Install necessary R packages
RUN R -e "install.packages('circlize')"
RUN R -e "install.packages('optparse')"
RUN R -e "install.packages('qpdf')"
RUN R -e "install.packages('pdftools')"

CMD ["R"]
