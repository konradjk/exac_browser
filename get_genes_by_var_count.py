import pymongo
import sys

num_genes = int(sys.argv[1])
db = pymongo.Connection().exac
tuples = []
for gene in db.genes.find():
    gene_id = gene['gene_id']
    num_variants = db.variants.find({'genes': gene_id}).count()
    tuples.append((gene_id, num_variants))
tuples = sorted(tuples, key=lambda x: x[1], reverse=True)
for gene_id, num_variants in tuples[:num_genes]:
    print gene_id