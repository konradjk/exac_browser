FROM gcr.io/exac-gnomad/exacbase

MAINTAINER MacArthur Lab

# python manage.py load_db \
# python manage.py load_db \
# python manage.py precalculate_metrics \
# python manage.py create_cache

RUN mkdir /var/exac_data;

# Fuse storage bucket
# ENTRYPOINT ["gcsfuse", "--implicit-dirs", "--key-file", "/var/www/docker/exac-gnomad-98cf40989f01.json", "exac", "/var/exac_data"]

# ENTRYPOINT ["./docker/mount.sh"]

CMD gcsfuse \
  --implicit-dirs \
  --key-file=/var/www/deploy/keys/exac-gnomad-30ea80400948.json \
  exac /var/exac_data && \
  ls /var/exac_data && \
  python manage.py load_db
# CMD ["/var/exac_data"]


# ENTRYPOINT ["python", "manage.py"]
# CMD ["load_gene_models"]
# CMD ["load_db"]