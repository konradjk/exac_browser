import json
import os
import pymongo
import gzip
from parsing import get_variants_from_sites_vcf, get_canonical_transcripts, \
    get_genes_from_gencode_gtf, get_transcripts_from_gencode_gtf, get_exons_from_gencode_gtf, \
    get_base_coverage_from_file, get_omim_associations, get_dbnsfp_info
import lookups
import xbrowse
from utils import *

from flask import Flask, request, session, g, redirect, url_for, abort, render_template, flash, jsonify
from flask.ext.compress import Compress
from flask_errormail import mail_on_500

from flask import Response
from collections import defaultdict
from werkzeug.contrib.cache import SimpleCache

from multiprocessing import Process
import glob
import time

ADMINISTRATORS = (
    'exac.browser.errors@gmail.com',
)

app = Flask(__name__)
mail_on_500(app, ADMINISTRATORS)
Compress(app)
app.config['COMPRESS_DEBUG'] = True
cache = SimpleCache()

EXAC_FILES_DIRECTORY = '../exac_test_data/'
REGION_LIMIT = 1E5
EXON_PADDING = 50
# Load default config and override config from an environment variable
app.config.update(dict(
    DB_HOST='localhost',
    DB_PORT=27017, 
    DB_NAME='exac', 
    DEBUG=True,
    SECRET_KEY='development key',

    SITES_VCFS=glob.glob(os.path.join(os.path.dirname(__file__), EXAC_FILES_DIRECTORY, 'sites*')),
    GENCODE_GTF=os.path.join(os.path.dirname(__file__), EXAC_FILES_DIRECTORY, 'gencode.gtf.gz'),
    CANONICAL_TRANSCRIPT_FILE=os.path.join(os.path.dirname(__file__), EXAC_FILES_DIRECTORY, 'canonical_transcripts.txt.gz'),
    OMIM_FILE=os.path.join(os.path.dirname(__file__), EXAC_FILES_DIRECTORY, 'omim_info.txt.gz'),
    BASE_COVERAGE_FILES=glob.glob(os.path.join(os.path.dirname(__file__), EXAC_FILES_DIRECTORY, 'coverage_split*')),
    DBNSFP_FILE=os.path.join(os.path.dirname(__file__), EXAC_FILES_DIRECTORY, 'dbNSFP2.6_gene.gz'),
))
GENE_CACHE_DIR = os.path.join(os.path.dirname(__file__), 'gene_cache')
GENES_TO_CACHE = {l.strip('\n') for l in open(os.path.join(os.path.dirname(__file__), 'genes_to_cache.txt'))}

def connect_db():
    """
    Connects to the specific database.
    """
    client = pymongo.MongoClient(host=app.config['DB_HOST'], port=app.config['DB_PORT'])
    return client[app.config['DB_NAME']]


def load_variants(sites_file, db, start_time):
    batch_size = 1000
    sites_vcf = gzip.open(sites_file)
    #size = os.path.getsize(sites_file)
    #progress = xbrowse.utils.get_progressbar(size, 'Loading Variants')
    current_entry = 0
    variants = []
    for variant in get_variants_from_sites_vcf(sites_vcf):
        current_entry += 1
        #progress.update(sites_vcf.fileobj.tell())
        variants.append(variant)
        if not current_entry % batch_size:
            db.variants.insert(variants, w=0)
            variants = []
            if not current_entry % 10000:
                print '%s up to %s (%s seconds so far)' % (sites_file, current_entry, (time.time() - start_time))
    if len(variants) > 0: db.variants.insert(variants, w=0)
    #progress.finish()


def load_coverage(coverage_fname, db, start_time):
    batch_size = 1000
    #size = os.path.getsize(coverage_file)
    coverage_file = gzip.open(coverage_fname)
    #progress = xbrowse.utils.get_progressbar(size, 'Parsing coverage')
    current_entry = 0
    bases = []
    for base_coverage in get_base_coverage_from_file(coverage_file):
        current_entry += 1
        #progress.update(coverage_file.fileobj.tell())
        bases.append(base_coverage)
        if not current_entry % batch_size:
            db.base_coverage.insert(bases, w=0)
            bases = []
            if not current_entry % 100000:
                print '%s up to %s (%s seconds so far)' % (coverage_fname, current_entry, (time.time() - start_time))
    if len(bases) > 0: db.base_coverage.insert(bases, w=0)
    #progress.finish()


def load_db():
    """
    Load the database
    """
    db = get_db()

    # Initialize database 
    # Don't need to explicitly create tables with mongo, just indices

    db.variants.drop()
    db.genes.drop()
    db.transcripts.drop()
    db.exons.drop()
    db.base_coverage.drop()

    # load coverage first; variant info will depend on coverage
    start_time = time.time()
    procs = []
    for fname in app.config['BASE_COVERAGE_FILES']:
        p = Process(target=load_coverage, args=(fname, db, start_time,))
        p.start()
        procs.append(p)
    [p.join() for p in procs]
    print 'Done loading. Took %s seconds' % (time.time() - start_time)

    db.base_coverage.ensure_index('xpos')

    # grab variants from sites VCF
    start_time = time.time()
    procs = []
    for fname in app.config['SITES_VCFS']:
        p = Process(target=load_variants, args=(fname, db, start_time,))
        p.start()
        procs.append(p)
    [p.join() for p in procs]
    print 'Done loading. Took %s seconds' % (time.time() - start_time)

    db.variants.ensure_index('xpos')
    db.variants.ensure_index('xstart')
    db.variants.ensure_index('xstop')
    db.variants.ensure_index('rsid')
    db.variants.ensure_index('genes')
    db.variants.ensure_index('transcripts')

    # grab genes from GTF
    gtf_file = gzip.open(app.config['GENCODE_GTF'])
    size = os.path.getsize(app.config['GENCODE_GTF'])
    #progress = xbrowse.utils.get_progressbar(size, 'Loading Genes')
    for gene in get_genes_from_gencode_gtf(gtf_file):
        db.genes.insert(gene, w=0)
        #progress.update(gtf_file.fileobj.tell())
    #progress.finish()
    gtf_file.close()

    db.genes.ensure_index('gene_id')
    db.genes.ensure_index('gene_name')
    db.genes.ensure_index('other_names')

    # and now transcripts
    gtf_file = gzip.open(app.config['GENCODE_GTF'])
    size = os.path.getsize(app.config['GENCODE_GTF'])
    #progress = xbrowse.utils.get_progressbar(size, 'Loading Transcripts')
    for transcript in get_transcripts_from_gencode_gtf(gtf_file):
        db.transcripts.insert(transcript, w=0)
        #progress.update(gtf_file.fileobj.tell())
    gtf_file.close()
    #progress.finish()

    db.transcripts.ensure_index('transcript_id')
    db.transcripts.ensure_index('gene_id')

    # Building up gene definitions
    gtf_file = gzip.open(app.config['GENCODE_GTF'])
    size = os.path.getsize(app.config['GENCODE_GTF'])
    #progress = xbrowse.utils.get_progressbar(size, 'Loading Exons')
    current_entry = 0
    exons = []
    for exon in get_exons_from_gencode_gtf(gtf_file):
        current_entry += 1
        exons.append(exon)
        db.exons.insert(exon, w=0)
        #progress.update(gtf_file.fileobj.tell())
        if not current_entry % 1000:
            db.base_coverage.insert(exons, w=0)
            exons = []
    gtf_file.close()
    #progress.finish()

    db.exons.ensure_index('exon_id')
    db.exons.ensure_index('transcript_id')
    db.exons.ensure_index('gene_id')

    canonical_transcript_file = gzip.open(app.config['CANONICAL_TRANSCRIPT_FILE'])
    size = os.path.getsize(app.config['CANONICAL_TRANSCRIPT_FILE'])
    #progress = xbrowse.utils.get_progressbar(size, 'Loading Canonical Transcripts')
    for gene, transcript in get_canonical_transcripts(canonical_transcript_file):
        gene = db.genes.find_one({
            'gene_id': gene
        })
        if not gene:
            continue
        gene['canonical_transcript'] = transcript
        db.genes.save(gene)
        #progress.update(canonical_transcript_file.fileobj.tell())
    canonical_transcript_file.close()
    #progress.finish()

    omim_file = gzip.open(app.config['OMIM_FILE'])
    size = os.path.getsize(app.config['OMIM_FILE'])
    #progress = xbrowse.utils.get_progressbar(size, 'Loading OMIM accessions')
    for fields in get_omim_associations(omim_file):
        if fields is None:
            continue
        gene, transcript, accession, description = fields
        gene = db.genes.find_one({
            'gene_id': gene
        })
        if not gene:
            continue
        gene['omim_accession'] = accession
        gene['omim_description'] = description
        db.genes.save(gene)
        #progress.update(omim_file.fileobj.tell())
    omim_file.close()
    #progress.finish()

    dbnsfp_file = gzip.open(app.config['DBNSFP_FILE'])
    size = os.path.getsize(app.config['DBNSFP_FILE'])
    #progress = xbrowse.utils.get_progressbar(size, 'Loading dbNSFP info')
    for dbnsfp_gene in get_dbnsfp_info(dbnsfp_file):
        gene = db.genes.find_one({
            'gene_id': dbnsfp_gene['ensembl_gene']
        })
        if not gene:
            continue
        gene['full_gene_name'] = dbnsfp_gene['gene_full_name']
        gene['other_names'] = []
        for other_name in dbnsfp_gene['gene_other_names']:
            gene['other_names'].append(other_name)
        db.genes.save(gene)
        #progress.update(dbnsfp_file.fileobj.tell())
    dbnsfp_file.close()
    #progress.finish()


def create_cache():
    """
    This is essentially a compile step that generates all cached resources.
    Creates files like autocomplete_entries.txt
    Should be run on every redeploy.
    """
    # create autocomplete_entries.txt
    autocomplete_strings = []
    for gene in get_db().genes.find():
        autocomplete_strings.append(gene['gene_name'])
        if 'other_names' in gene:
            autocomplete_strings.extend(gene['other_names'])
    f = open(os.path.join(os.path.dirname(__file__), 'autocomplete_strings.txt'), 'w')
    for s in autocomplete_strings:
        f.write(s+'\n')
    f.close()

    # create static gene pages for genes in
    if not os.path.exists(GENE_CACHE_DIR):
        os.makedirs(GENE_CACHE_DIR)

    # get list of genes ordered by num_variants
    for gene_id in GENES_TO_CACHE:
        try:
            page_content = get_gene_page_content(gene_id)
        except Exception:
            print Exception
            continue
        f = open(os.path.join(GENE_CACHE_DIR, '{}.html'.format(gene_id)), 'w')
        f.write(page_content)
        f.close()


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
    return render_template('homepage.html')


@app.route('/autocomplete/<query>')
def awesome_autocomplete(query):
    if not hasattr(g, 'autocomplete_strings'):
        g.autocomplete_strings = [s.strip() for s in open(os.path.join(os.path.dirname(__file__), 'autocomplete_strings.txt'))]
    suggestions = lookups.get_awesomebar_suggestions(g, query)
    return Response(json.dumps([{'value': s} for s in suggestions]),  mimetype='application/json')


@app.route('/awesome')
def awesome():
    db = get_db()
    query = request.args.get('query')
    datatype, identifier = lookups.get_awesomebar_result(db, query)

    print "Searched for %s %s" % (datatype, identifier)
    if datatype == 'gene':
        return redirect('/gene/{}'.format(identifier))
    elif datatype == 'transcript':
        return redirect('/transcript/{}'.format(identifier))
    elif datatype == 'variant':
        return redirect('/variant/{}'.format(identifier))
    elif datatype == 'region':
        return redirect('/region/{}'.format(identifier))
    elif datatype == 'dbsnp_variant_set':
        return redirect('/dbsnp/{}'.format(identifier))
    elif datatype == 'error':
        return redirect('/error/{}'.format(identifier))
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

    # pos, ref, alt = get_minimal_representation(pos, ref, alt)
    xpos = xbrowse.get_xpos(chrom, pos)
    variant = lookups.get_variant(db, xpos, ref, alt)

    if variant is None:
        variant = {
            'chrom': chrom,
            'pos': pos,
            'xpos': xpos,
            'ref': ref,
            'alt': alt
        }
    consequences = None
    ordered_csqs = None
    if 'vep_annotations' in variant:
        variant['vep_annotations'] = order_vep_by_csq(variant['vep_annotations'])  # Adds major_consequence
        ordered_csqs = [x['major_consequence'] for x in variant['vep_annotations']]
        ordered_csqs = reduce(lambda x, y: ','.join([x, y]) if y not in x else x, ordered_csqs).split(',') # Close but not quite there
        consequences = defaultdict(lambda: defaultdict(list))
        for annotation in variant['vep_annotations']:
            consequences[annotation['major_consequence']][annotation['Gene']].append(annotation)
    base_coverage = lookups.get_coverage_for_bases(db, xpos, xpos + len(ref) - 1)
    any_covered = any([x['has_coverage'] for x in base_coverage])
    return render_template(
        'variant.html',
        variant=variant,
        base_coverage=base_coverage,
        consequences=consequences,
        any_covered=any_covered,
        ordered_csqs=ordered_csqs
    )


@app.route('/gene/<gene_id>')
def gene_page(gene_id):
    if gene_id in GENES_TO_CACHE:
        return open(os.path.join(GENE_CACHE_DIR, '{}.html'.format(gene_id))).read()
    else:
        return get_gene_page_content(gene_id)

def get_gene_page_content(gene_id):
    db = get_db()
    gene = lookups.get_gene(db, gene_id)
    if gene is None:
        abort(404)
    cache_key = 't-gene-{}'.format(gene_id)
    t = cache.get(cache_key)
    if t is None:
        variants_in_gene = lookups.get_variants_in_gene(db, gene_id)
        transcripts_in_gene = lookups.get_transcripts_in_gene(db, gene_id)

        # Get some canonical transcript and corresponding info
        transcript_id = gene['canonical_transcript']
        transcript = lookups.get_transcript(db, transcript_id)
        variants_in_transcript = lookups.get_variants_in_transcript(db, transcript_id)
        coverage_stats = lookups.get_coverage_for_transcript(db, transcript['xstart'] - EXON_PADDING, transcript['xstop'] + EXON_PADDING)
        add_transcript_coordinate_to_variants(db, variants_in_transcript, transcript_id)

        t = render_template(
            'gene.html',
            gene=gene,
            transcript=transcript,
            variants_in_gene=variants_in_gene,
            variants_in_transcript=variants_in_transcript,
            transcripts_in_gene=transcripts_in_gene,
            coverage_stats=coverage_stats
        )
        cache.set(cache_key, t, timeout=1000*60)
    return t


@app.route('/transcript/<transcript_id>')
def transcript_page(transcript_id):
    db = get_db()
    transcript = lookups.get_transcript(db, transcript_id)

    cache_key = 't-transcript-{}'.format(transcript_id)
    t = cache.get(cache_key)
    if t is None: 
    
        gene = lookups.get_gene(db, transcript['gene_id'])
        gene['transcripts'] = lookups.get_transcripts_in_gene(db, transcript['gene_id'])
        variants_in_transcript = lookups.get_variants_in_transcript(db, transcript_id)

        coverage_stats = lookups.get_coverage_for_transcript(db, transcript['xstart'] - EXON_PADDING, transcript['xstop'] + EXON_PADDING)

        add_transcript_coordinate_to_variants(db, variants_in_transcript, transcript_id)

        t = render_template(
            'transcript.html',
            transcript=transcript,
            transcript_json=json.dumps(transcript),
            variants_in_transcript=variants_in_transcript,
            variants_in_transcript_json=json.dumps(variants_in_transcript),
            coverage_stats=coverage_stats,
            coverage_stats_json=json.dumps(coverage_stats),
            gene=gene,
            gene_json=json.dumps(gene),
        )
        cache.set(cache_key, t, timeout=1000*60)
    return t


@app.route('/region/<region_id>')
def region_page(region_id):
    db = get_db()
    region = region_id.split('-')
    cache_key = 't-region-{}'.format(region_id)
    t = cache.get(cache_key)
    if t is None:
        chrom = region[0]
        start = None
        stop = None
        if len(region) == 3:
            chrom, start, stop = region
            start = int(start)
            stop = int(stop)
        if start is None or stop - start > REGION_LIMIT:
            return render_template(
                'region.html',
                genes_in_region=None,
                variants_in_region=None,
                chrom=chrom,
                start=start,
                stop=stop,
                coverage=None
            )
        genes_in_region = lookups.get_genes_in_region(db, chrom, start, stop)
        variants_in_region = lookups.get_variants_in_region(db, chrom, start, stop)
        xstart = xbrowse.get_xpos(chrom, start)
        xstop = xbrowse.get_xpos(chrom, stop)
        coverage_array = lookups.get_coverage_for_bases(db, xstart, xstop)
        t = render_template(
            'region.html',
            genes_in_region=genes_in_region,
            variants_in_region=variants_in_region,
            chrom=chrom,
            start=start,
            stop=stop,
            coverage=coverage_array
        )
    return t


@app.route('/dbsnp/<rsid>')
def dbsnp_page(rsid):
    db = get_db()
    variants = lookups.get_variants_by_rsid(db, rsid)
    chrom = None
    start = None
    stop = None
    return render_template(
        'region.html',
        rsid=rsid,
        variants_in_region=variants,
        chrom=chrom,
        start=start,
        stop=stop,
        coverage=None,
        genes_in_region=None
    )


@app.route('/error/<query>')
@app.errorhandler(404)
def error_page(query):
    unsupported = "TTN" if query in lookups.UNSUPPORTED_QUERIES else None

    return render_template(
        'error.html',
        query=query,
        unsupported=unsupported
    )


@app.route('/downloads')
def downloads_page():
    return render_template('downloads.html')


@app.route('/about')
def about_page():
    return render_template('about.html')


@app.route('/participants')
def participants_page():
    return render_template('about.html')


@app.route('/terms')
def terms_page():
    return render_template('terms.html')


if __name__ == "__main__":
    app.run()
