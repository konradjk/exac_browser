FROM gcr.io/exac-gnomad/exacbase

MAINTAINER MacArthur Lab

COPY . /var/www/

RUN mkdir /var/exac_data;

EXPOSE 80

CMD gcsfuse \
  --implicit-dirs \
  --key-file=/var/www/deploy/keys/exac-gnomad-30ea80400948.json \
  exac /var/exac_data && \
  ls /var/exac_data && \
  python exac.py -h 0.0.0.0 -p 80
