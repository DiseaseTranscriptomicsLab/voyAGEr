FROM r-base:4.1.2
MAINTAINER Disease Transcriptomics Lab <imm-nmorais@medicina.ulisboa.pt>

RUN sed -i 's|http://deb.debian.org|http://ftp.us.debian.org|g' /etc/apt/sources.list && \
    apt-get update && apt-get -y upgrade && apt-get -y autoremove
RUN apt-get install -y curl
RUN apt-get install -y libcurl4-openssl-dev
RUN apt-get install -y libssl-dev
RUN R -e 'install.packages("shiny")'
RUN R -e 'install.packages("shinycssloaders")'
RUN R -e 'install.packages("ggplot2")'
RUN R -e 'install.packages("fontawesome")'
RUN R -e 'install.packages("manipulateWidget")'
RUN R -e 'install.packages("shinycssloaders")'
RUN R -e 'install.packages("viridis")'
RUN R -e 'install.packages("RSQLite")'
RUN R -e 'install.packages("manipulateWidget")'
RUN R -e 'install.packages("highcharter")'
RUN R -e 'install.packages("tm")'
RUN R -e 'install.packages("reshape2")'
RUN R -e 'install.packages("plyr")'
RUN R -e 'install.packages("htmltools")'
RUN R -e 'install.packages("RColorBrewer")'
RUN R -e 'install.packages("shinyBS")'
RUN R -e 'install.packages("shinythemes")'
RUN R -e 'install.packages("shinyWidgets")'
RUN R -e 'install.packages("reshape2")'
RUN apt-get install -y libxml2-dev
RUN R -e 'install.packages("highcharter")'
RUN R -e 'install.packages("DT")'
RUN R -e 'install.packages("tm")'
RUN R -e 'install.packages("circlize")'

RUN R -e "install.packages(c('withr'), repos='https://cloud.r-project.org/')"
RUN R -e "withr::with_makevars(c(PKG_CPPFLAGS='-DHTTP_MAX_HEADER_SIZE=0x7fffffff'), {install.packages(c('shiny'), repos='https://cloud.r-project.org/')}, assignment = '+=')"

WORKDIR /home/app

# copy the app directory into the image
COPY . /home/app

#COPY ./Rprofile.site /usr/lib/R/etc/

#EXPOSE 3838
CMD ["R", "-e", "shiny::runApp(host='0.0.0.0', port=3838)"]
