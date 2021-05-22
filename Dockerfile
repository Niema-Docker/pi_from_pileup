# Minimal Docker image for pi_from_pileup using Alpine base
FROM alpine:latest
MAINTAINER Niema Moshiri <niemamoshiri@gmail.com>

# install pi_from_pileup
RUN apk update && \
    apk add bash g++ && \
    wget "https://raw.githubusercontent.com/Niema-Docker/pi_from_pileup/main/pi_from_pileup.cpp" && \
    g++ -O3 -o /usr/local/bin/pi_from_pileup pi_from_pileup.cpp && \
    rm pi_from_pileup.cpp
