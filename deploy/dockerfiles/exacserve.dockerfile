FROM gcr.io/exac-gnomad/exacbase

MAINTAINER MacArthur Lab

EXPOSE 80

ENTRYPOINT ["python"]
CMD ["exac.py", "-h 0.0.0.0", "-p 80"]

