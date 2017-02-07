FROM gcr.io/exac-gnomad/exacbase

MAINTAINER MacArthur Lab

CMD pip install -U crcmod && \
  gsutil -m cp -r gs://exac/combined_bams_v3 /readviz/ && \
  gsutil cp gs://exac/gencode.v19.sorted.bed /readviz/ && \
  gsutil cp gs://exac/gencode.v19.sorted.bed.idx /readviz/
