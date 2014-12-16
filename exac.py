import itertools
import json
import os
import pymongo
import pysam
import gzip
from parsing import get_variants_from_sites_vcf, get_canonical_transcripts, \
    get_genes_from_gencode_gtf, get_transcripts_from_gencode_gtf, get_exons_from_gencode_gtf, \
    get_base_coverage_from_file, get_omim_associations, get_dbnsfp_info, get_snp_from_dbsnp_file
import lookups
import xbrowse
import random
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
    LOAD_DB_PARALLEL_PROCESSES = 4,  # contigs assigned to threads, so good to make this a factor of 24 (eg. 2,3,4,6,8)
    #SITES_VCFS=glob.glob(os.path.join(os.path.dirname(__file__), EXAC_FILES_DIRECTORY, 'sites*')),
    SITES_VCFS=glob.glob(os.path.join(os.path.dirname(__file__), EXAC_FILES_DIRECTORY, 'ExAC.r0.2.*.vep.vcf.gz')),
    GENCODE_GTF=os.path.join(os.path.dirname(__file__), EXAC_FILES_DIRECTORY, 'gencode.gtf.gz'),
    CANONICAL_TRANSCRIPT_FILE=os.path.join(os.path.dirname(__file__), EXAC_FILES_DIRECTORY, 'canonical_transcripts.txt.gz'),
    OMIM_FILE=os.path.join(os.path.dirname(__file__), EXAC_FILES_DIRECTORY, 'omim_info.txt.gz'),
    BASE_COVERAGE_FILES=glob.glob(os.path.join(os.path.dirname(__file__), EXAC_FILES_DIRECTORY, 'Panel.*.coverage.txt.gz')),
    DBNSFP_FILE=os.path.join(os.path.dirname(__file__), EXAC_FILES_DIRECTORY, 'dbNSFP2.6_gene.gz'),
    DBSNP_FILE=os.path.join(os.path.dirname(__file__), EXAC_FILES_DIRECTORY, 'snp141.txt.gz') # wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/snp141.txt.gz
))
GENE_CACHE_DIR = os.path.join(os.path.dirname(__file__), 'gene_cache')
GENES_TO_CACHE = {l.strip('\n') for l in open(os.path.join(os.path.dirname(__file__), 'genes_to_cache.txt'))}


def connect_db():
    """
    Connects to the specific database.
    """
    client = pymongo.MongoClient(host=app.config['DB_HOST'], port=app.config['DB_PORT'])
    return client[app.config['DB_NAME']]



def load_base_coverage():
    def coverage_parser(coverage_files, start_time):
        for coverage_fname in coverage_files:
            with gzip.open(coverage_fname) as coverage_file:
                #progress = xbrowse.utils.get_progressbar(size, 'Parsing coverage')
                #current_entry = 0
                for base_coverage in get_base_coverage_from_file(coverage_file):
                    #current_entry += 1
                    yield base_coverage
                    #progress.update(coverage_file.fileobj.tell())
                    #if current_entry % 1000000 == 0:
                    #  print '%s up to %s (%s seconds so far)' % (coverage_fname, current_entry, int(time.time()-start_time))
        #progress.finish()

    def load_coverage(coverage_files, db, start_time):
        db.base_coverage.insert(coverage_parser(coverage_files, start_time), w=0)

    db = get_db()
    db.base_coverage.drop()
    print("Dropped db.base_coverage")
    # load coverage first; variant info will depend on coverage
    start_time = time.time()
    db.base_coverage.ensure_index('xpos')

    procs = []
    num_procs = app.config['LOAD_DB_PARALLEL_PROCESSES']
    random.shuffle(app.config['BASE_COVERAGE_FILES'])
    for i in range(num_procs):
        coverage_files_subset = app.config['BASE_COVERAGE_FILES'][i::num_procs]
        p = Process(target=load_coverage, args=(coverage_files_subset, db, start_time,))
        p.start()
        procs.append(p)
        print("Created process %(i)d to load files %(coverage_files_subset)s from %(coverage_files_subset)s" % locals())
    [p.join() for p in procs]

    print 'Done loading coverage. Took %s seconds' % (time.time() - start_time)


def load_variants_file():

    def variant_parser(tabix_file, contigs_to_load, start_time):
        for contig in contigs_to_load:
            vcf_header_iterator = tabix_file.header
            vcf_records_iterator = tabix_file.fetch(contig, 0, 10**9, multiple_iterators=True)

            current_entry = 0
            for variant in get_variants_from_sites_vcf(itertools.chain(vcf_header_iterator, vcf_records_iterator)):
                current_entry += 1
                yield variant

                #if current_entry % 10000 == 0:
                #    print 'Loading %s, at chrom %s: %s  (%s seconds so far)' % (
                #        os.path.basename(tabix_file.filename), variant['chrom'], variant['pos'], int(time.time()-start_time))
            print("Finished loading %s variants from chrom %s" % (current_entry, contig))


    # Passes a generator to db.variants.insert(..), and lets it decide how many variants to insert at a time
    # Based on http://api.mongodb.org/python/current/examples/bulk.html
    def insert_variants(tabix_file, db, contigs_to_load, start_time):
        db.variants.insert(variant_parser(tabix_file, contigs_to_load, start_time), w=0)

    db = get_db()
    db.variants.drop()
    print("Dropped db.variants")

    # grab variants from sites VCF
    start_time = time.time()
    db.variants.ensure_index('xpos')
    db.variants.ensure_index('xstart')
    db.variants.ensure_index('xstop')
    db.variants.ensure_index('rsid')
    db.variants.ensure_index('genes')
    db.variants.ensure_index('transcripts')

    for fname in app.config['SITES_VCFS']:
        procs = []
        f = pysam.TabixFile(fname)
        random.shuffle(f.contigs)
        num_procs = app.config['LOAD_DB_PARALLEL_PROCESSES']
        for i in range(num_procs):
            contigs_for_this_proc = f.contigs[i::num_procs]
            p = Process(target=insert_variants, args=(f, db, contigs_for_this_proc, start_time))
            p.start()
            procs.append(p)
            print("Created process %(i)d to load contigs %(contigs_for_this_proc)s from %(fname)s " % locals())
        [p.join() for p in procs]
        f.close()
    print 'Done loading variants. Took %s seconds' % (time.time() - start_time)


def load_gene_models():
    db = get_db()

    db.genes.drop()
    db.transcripts.drop()
    db.exons.drop()

    start_time = time.time()

    #size = os.path.getsize(app.config['CANONICAL_TRANSCRIPT_FILE'])
    canonical_transcripts = {}
    with gzip.open(app.config['CANONICAL_TRANSCRIPT_FILE']) as canonical_transcript_file:
        #progress = xbrowse.utils.get_progressbar(size, 'Loading Canonical Transcripts')
        for gene, transcript in get_canonical_transcripts(canonical_transcript_file):
            canonical_transcripts[gene] = transcript
            #progress.update(canonical_transcript_file.fileobj.tell())
        #progress.finish()

    #size = os.path.getsize(app.config['OMIM_FILE'])
    omim_annotations = {}
    with gzip.open(app.config['OMIM_FILE']) as omim_file:
        #progress = xbrowse.utils.get_progressbar(size, 'Loading OMIM accessions')
        for fields in get_omim_associations(omim_file):
            if fields is None:
                continue
            gene, transcript, accession, description = fields
            omim_annotations[gene] = (accession, description)
            #progress.update(omim_file.fileobj.tell())
        #progress.finish()

    #size = os.path.getsize(app.config['DBNSFP_FILE'])
    dbnsfp_info = {}
    with gzip.open(app.config['DBNSFP_FILE']) as dbnsfp_file:
        #progress = xbrowse.utils.get_progressbar(size, 'Loading dbNSFP info')
        for dbnsfp_gene in get_dbnsfp_info(dbnsfp_file):
            other_names = [other_name.upper() for other_name in dbnsfp_gene['gene_other_names']]
            dbnsfp_info[dbnsfp_gene['ensembl_gene']] = (dbnsfp_gene['gene_full_name'], other_names)
            #progress.update(dbnsfp_file.fileobj.tell())
        #progress.finish()

    print 'Done loading metadata. Took %s seconds' % (time.time() - start_time)

    # grab genes from GTF
    start_time = time.time()

    #size = os.path.getsize(app.config['GENCODE_GTF'])
    #progress = xbrowse.utils.get_progressbar(size, 'Loading Genes')
    with gzip.open(app.config['GENCODE_GTF']) as gtf_file:
        for gene in get_genes_from_gencode_gtf(gtf_file):
            gene_id = gene['gene_id']
            if gene_id in canonical_transcripts:
                gene['canonical_transcript'] = canonical_transcripts[gene_id]
            if gene_id in omim_annotations:
                gene['omim_accession'] = omim_annotations[gene_id][0]
                gene['omim_description'] = omim_annotations[gene_id][1]
            if gene_id in dbnsfp_info:
                gene['full_gene_name'] = dbnsfp_info[gene_id][0]
                gene['other_names'] = dbnsfp_info[gene_id][1]
            db.genes.insert(gene, w=0)
            #progress.update(gtf_file.fileobj.tell())
        #progress.finish()

    print 'Done loading genes. Took %s seconds' % (time.time() - start_time)

    start_time = time.time()
    db.genes.ensure_index('gene_id')
    db.genes.ensure_index('gene_name_upper')
    db.genes.ensure_index('gene_name')
    db.genes.ensure_index('other_names')
    db.genes.ensure_index('xstart')
    db.genes.ensure_index('xstop')
    print 'Done indexing gene table. Took %s seconds' % (time.time() - start_time)

    # and now transcripts
    start_time = time.time()
    #size = os.path.getsize(app.config['GENCODE_GTF'])
    #progress = xbrowse.utils.get_progressbar(size, 'Loading Transcripts')
    with gzip.open(app.config['GENCODE_GTF']) as gtf_file:
        db.transcripts.insert((transcript for transcript in get_transcripts_from_gencode_gtf(gtf_file)), w=0)
    #progress.update(gtf_file.fileobj.tell())
    #progress.finish()
    print 'Done loading transcripts. Took %s seconds' % (time.time() - start_time)

    start_time = time.time()
    db.transcripts.ensure_index('transcript_id')
    db.transcripts.ensure_index('gene_id')
    print 'Done indexing transcript table. Took %s seconds' % (time.time() - start_time)

    # Building up gene definitions
    start_time = time.time()
    gtf_file = gzip.open(app.config['GENCODE_GTF'])
    size = os.path.getsize(app.config['GENCODE_GTF'])
    #progress = xbrowse.utils.get_progressbar(size, 'Loading Exons')
    current_entry = 0
    exons = []
    for exon in get_exons_from_gencode_gtf(gtf_file):
        current_entry += 1
        exons.append(exon)
        #progress.update(gtf_file.fileobj.tell())
        if not current_entry % 1000:
            db.exons.insert(exons, w=0)
            exons = []
    if len(exons) > 0: db.exons.insert(exons, w=0)
    gtf_file.close()
    #progress.finish()
    print 'Done loading exons. Took %s seconds' % (time.time() - start_time)

    start_time = time.time()
    db.exons.ensure_index('exon_id')
    db.exons.ensure_index('transcript_id')
    db.exons.ensure_index('gene_id')
    print 'Done indexing exon table. Took %s seconds' % (time.time() - start_time)


def load_dbsnp():
    db = get_db()

    db.dbsnp.drop()

    print "Loading dbsnp from %s " % app.config['DBSNP_FILE']
    start_time = time.time()
    db.dbsnp.ensure_index('rsid')
    with gzip.open(app.config['DBSNP_FILE']) as dbsnp_file:
        db.dbsnp.insert((snp for snp in get_snp_from_dbsnp_file(dbsnp_file)), w=0)
    print 'Done loading dbSNP. Took %s seconds' % (time.time() - start_time)

    #start_time = time.time()
    #db.dbsnp.ensure_index('rsid')
    #print 'Done indexing dbSNP table. Took %s seconds' % (time.time() - start_time)


def load_db():
    """
    Load the database
    """
    # Initialize database
    # Don't need to explicitly create tables with mongo, just indices
    procs = []
    for load_function in [load_dbsnp, load_variants_file, load_base_coverage, load_gene_models]:
        p = Process(target=load_function)
        p.start()
        procs.append(p)
        print("Started process for: " + load_function.__name__)

    [p.join() for p in procs]


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
    for s in sorted(autocomplete_strings):
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

    print "Searched for %s: %s" % (datatype, identifier)
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
            ordered_csqs = reduce(lambda x, y: ','.join([x, y]) if y not in x else x, ordered_csqs, '').split(',') # Close but not quite there
            consequences = defaultdict(lambda: defaultdict(list))
            for annotation in variant['vep_annotations']:
                annotation['HGVS'] = get_proper_hgvs(annotation)
                consequences[annotation['major_consequence']][annotation['Gene']].append(annotation)
        base_coverage = lookups.get_coverage_for_bases(db, xpos, xpos + len(ref) - 1)
        any_covered = any([x['has_coverage'] for x in base_coverage])
        print 'Rendering variant: %s' % variant_str
        return render_template(
            'variant.html',
            variant=variant,
            base_coverage=base_coverage,
            consequences=consequences,
            any_covered=any_covered,
            ordered_csqs=ordered_csqs
        )
    except Exception, e:
        print 'Failed on variant:', variant_str, '; Error=', e
        abort(404)


@app.route('/gene/<gene_id>')
def gene_page(gene_id):
    if gene_id in GENES_TO_CACHE:
        return open(os.path.join(GENE_CACHE_DIR, '{}.html'.format(gene_id))).read()
    else:
        return get_gene_page_content(gene_id)


def get_gene_page_content(gene_id):
    db = get_db()
    try:
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
        print 'Rendering gene: %s' % gene_id
        return t
    except Exception, e:
        print 'Failed on gene:', gene_id, ';Error=', e
        abort(404)


@app.route('/transcript/<transcript_id>')
def transcript_page(transcript_id):
    db = get_db()
    try:
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
        print 'Rendering transcript: %s' % transcript_id
        return t
    except Exception, e:
        print 'Failed on transcript:', transcript_id, ';Error=', e
        abort(404)


@app.route('/region/<region_id>')
def region_page(region_id):
    db = get_db()
    try:
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
        print 'Rendering region: %s' % region_id
        return t
    except Exception, e:
        print 'Failed on region:', region_id, ';Error=', e
        abort(404)


@app.route('/dbsnp/<rsid>')
def dbsnp_page(rsid):
    db = get_db()
    try:
        variants = lookups.get_variants_by_rsid(db, rsid)
        chrom = None
        start = None
        stop = None
        print 'Rendering rsid: %s' % rsid
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
    except Exception, e:
        print 'Failed on rsid:', rsid, ';Error=', e
        abort(404)


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


@app.route('/contact')
def contact_page():
    return render_template('contact.html')


@app.route('/faq')
def faq_page():
    return render_template('faq.html')


if __name__ == "__main__":
    app.run()
