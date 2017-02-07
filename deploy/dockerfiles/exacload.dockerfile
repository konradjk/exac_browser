FROM gcr.io/exac-gnomad/exacbase

MAINTAINER MacArthur Lab

COPY . /var/www/

RUN mkdir /var/exac_data;

CMD gcsfuse \
  --implicit-dirs \
  --key-file=/var/www/deploy/keys/exac-gnomad-30ea80400948.json \
  exac /var/exac_data && \
  ls /var/exac_data && \
  python manage.py load_db
  # python manage.py load_base_coverage
  # python manage.py load_constraint_information
  # python manage.py load_constraint_information && \
  # python manage.py load_mnps && \
  # python manage.py create_cache && \
  # python manage.py precalculate_metrics
