import json
import os
import re
import pymongo
import gzip
from parsing import get_variants_from_sites_vcf, get_genes_from_gencode_gtf, get_transcripts_from_gencode_gtf
import lookups
import xbrowse
import copy
#from xbrowse.annotation.vep_annotations import get_vep_annotations_from_vcf

from flask import Flask, request, session, g, redirect, url_for, abort, render_template, flash, jsonify
from flask import Response


app = Flask(__name__)

# Load default config and override config from an environment variable
app.config.update(dict(
    DB_HOST='localhost',
    DB_PORT=27017, 
    DB_NAME='exac', 
    DEBUG = True,
    SECRET_KEY = 'development key',

    SITES_VCF = os.path.join(os.path.dirname(__file__), '../exac_anno.vep_0001.vep.vcf.gz'),
    GENCODE_GTF = os.path.join(os.path.dirname(__file__), '../gencode.v19.annotation.gtf.gz'),

))


def connect_db():
    """
    Connects to the specific database.
    """
    client = pymongo.MongoClient(host=app.config['DB_HOST'], port=app.config['DB_PORT'])
    return client[app.config['DB_NAME']]


def load_db():
    """
    Load the database
    """
    db = get_db()

    # Initialize database 
    # Don't need to explicitly create tables with mongo, just indices
    db.variants.remove()
    db.genes.remove()
    db.transcripts.remove()

    db.variants.ensure_index('xpos')
    db.variants.ensure_index('genes')
    db.variants.ensure_index('trancripts')

    db.genes.ensure_index('gene_id')
    db.genes.ensure_index('gene_name')

    db.transcripts.ensure_index('transcript_id')

    # grab variants from sites VCF
    sites_vcf = gzip.open(app.config['SITES_VCF'])
    size = os.path.getsize(app.config['SITES_VCF'])
    progress = xbrowse.utils.get_progressbar(size, 'Loading Variants')
    for variant in get_variants_from_sites_vcf(sites_vcf):
        db.variants.insert(variant)
        progress.update(sites_vcf.fileobj.tell())

    # grab genes from GTF
    gtf_file = gzip.open(app.config['GENCODE_GTF'])
    size = os.path.getsize(app.config['GENCODE_GTF'])
    progress = xbrowse.utils.get_progressbar(size, 'Loading Genes')
    for gene in get_genes_from_gencode_gtf(gtf_file):
        db.genes.insert(gene)
        progress.update(gtf_file.fileobj.tell())
    gtf_file.close()

    # and now transcripts
    gtf_file = gzip.open(app.config['GENCODE_GTF'])
    size = os.path.getsize(app.config['GENCODE_GTF'])
    progress = xbrowse.utils.get_progressbar(size, 'Loading Transcripts')
    for transcript in get_transcripts_from_gencode_gtf(gtf_file):
        db.transcripts.insert(transcript)
        progress.update(gtf_file.fileobj.tell())
    gtf_file.close()


def get_db():
    """
    Opens a new database connection if there is none yet for the
    current application context.
    """
    if not hasattr(g, 'db_conn'):
        g.db_conn = connect_db()
    return g.db_conn


# @app.teardown_appcontext
# def close_db(error):
#     """Closes the database again at the end of the request."""
#     if hasattr(g, 'db_conn'):
#         g.db_conn.close()


@app.route('/')
def homepage():
    db = get_db()
    return render_template('homepage.html')


@app.route('/autocomplete/<query>')
def awesome_autocomplete(query):
    db = get_db()
    suggestions = lookups.get_awesomebar_suggestions(db, query)
    return Response(json.dumps([{'value': s} for s in suggestions]),  mimetype='application/json')


@app.route('/awesome')
def awesome():
    db = get_db()
    query = request.args.get('query')
    datatype, identifier = lookups.get_awesomebar_result(db, query)
    if datatype == 'gene':
        return redirect('/gene/{}'.format(identifier))
    elif datatype == 'variant':
        return redirect('/variant/{}'.format(identifier))
    elif datatype == 'region':
        raise NotImplementedError
    else:
        raise Exception


@app.route('/variant/<variant_str>')
def variant_page(variant_str):
    db = get_db()
    try:
        chrom, pos, ref, alt = variant_str.split('-')
        pos = int(pos)
    except ValueError:
        abort(404)
    xpos = xbrowse.get_xpos(chrom, pos)
    variant = lookups.get_variant(db, xpos, ref, alt)
    return render_template('variant.html', variant=variant)


@app.route('/gene/<gene_id>')
def gene_page(gene_id):
    db = get_db()
    gene = lookups.get_gene(db, gene_id)
    variants_in_gene = lookups.get_variants_in_gene(db, gene_id)
    return render_template('gene.html', gene=gene, variants_in_gene=variants_in_gene)


@app.route('/transcript/<transcript_id>')
def transcript_page(transcript_id):
    db = get_db()
    transcript = lookups.get_transcript(db, transcript_id)
    variants_in_transcript = lookups.get_variants_in_transcript(db, transcript_id)
    return render_template('transcript.html', transcript=transcript, variants_in_transcript=variants_in_transcript)


@app.route('/howtouse')
def howtouse_page():
    return render_template('howtouse.html')


@app.route('/downloads')
def downloads_page():
    return render_template('downloads.html')


@app.route('/contact')
def contact_page():
    return render_template('contact.html')


if __name__ == "__main__":
    app.run()
