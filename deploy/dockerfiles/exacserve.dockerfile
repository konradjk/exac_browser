FROM gcr.io/exac-gnomad/exacbase

MAINTAINER MacArthur Lab

COPY . /var/www/

RUN mkdir /var/exac_data;

EXPOSE 80

ENTRYPOINT ["python"]
CMD ["exac.py", "-h 0.0.0.0", "-p 80"]

