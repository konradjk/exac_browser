FROM python:2.7.12

MAINTAINER MacArthur Lab

# ENV GNOMAD_DB_HOST=gnomaddb GNOMAD_DB_PORT=27017 GNOMAD_DATA=../exac_data/ READ_VIZ=/mongo/readviz

COPY . /var/www
WORKDIR /var/www

RUN pip install -r requirements.txt