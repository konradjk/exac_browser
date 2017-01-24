FROM macarthurlab/exacbase

MAINTAINER MacArthur Lab

# ENV EXAC_DATA=gs://exac/170122_exacv1_bundle/ READ_VIZ=/mongo/readviz
ENV EXAC_DATA=/Users/msolomon/Projects/exacg/exacv1_bundle/ READ_VIZ=/mongo/readviz

# python manage.py load_db \
# python manage.py load_db \
# python manage.py precalculate_metrics \
# python manage.py create_cache

ENTRYPOINT ["python", "manage.py"]
CMD ["load_gene_models"]
# CMD ["load_db"]