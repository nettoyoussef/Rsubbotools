FROM rocker/verse:4.4.2

LABEL Name=general Version=0.0.1

ENV LISTEN_PORT=5001
EXPOSE 5001

# copy Rprofile configuration with favorite CRAN repo to 
# home directory
COPY ./.Rprofile ~/.Rprofile 

RUN echo "$RENV_PATHS_CACHE"

RUN mkdir -p ~/.local/share/renv/

################################

# Install RcppGSL dependencies

# To install gsl
# and the gdb debugger (gcc)
RUN apt-get update && apt-get install -y --force-yes --allow-unauthenticated --no-install-recommends --no-upgrade \      
    libgsl-dev  \
    gdb

################################


################################

# Install subbotools 

#RUN wget -P /tmp/ ftp://cafed.sssup.it/packages/subbotools-1.3.0.tar.gz && \
#    cd /tmp && \
#    tar xvzf subbotools-1.3.0.tar.gz && \
#    cd subbotools-1.3.0 && \
#    ./configure && \
#    make && \
#    make install

# tive que usar uma versao editada do pacote ja que a versao disponivel no site
# de Santanna esta com problemas nos arquivos do laplafit e laplaafit (Data e Size
# precisam ser declarados com extern)
COPY subbotools-1.3.0.tar.gz /tmp/subbotools-1.3.0.tar.gz 

RUN cd /tmp/ && \
    tar xvzf ./subbotools-1.3.0.tar.gz && \
    cd subbotools-1.3.0 && \
    ./configure && \
    make && \
    make install



################################

# restore the lockfile for the container
RUN R -e 'install.packages("renv")'

CMD ["/bin/bash"]

