import json
import os
import pymongo
import gzip
from parsing import get_variants_from_sites_vcf, get_canonical_transcripts, \
    get_genes_from_gencode_gtf, get_transcripts_from_gencode_gtf, get_exons_from_gencode_gtf, \
    get_base_coverage_from_file
import lookups
import xbrowse
from utils import *

from flask import Flask, request, session, g, redirect, url_for, abort, render_template, flash, jsonify
from flask import Response
from collections import defaultdict

app = Flask(__name__)

EXAC_FILES_DIRECTORY = '../exac_test_data/'
REGION_LIMIT = 1E6
# Load default config and override config from an environment variable
app.config.update(dict(
    DB_HOST='localhost',
    DB_PORT=27017, 
    DB_NAME='exac', 
    DEBUG=True,
    SECRET_KEY='development key',

    SITES_VCF=os.path.join(os.path.dirname(__file__), EXAC_FILES_DIRECTORY, 'sites_file.vcf.gz'),
    FULL_VCF=os.path.join(os.path.dirname(__file__), EXAC_FILES_DIRECTORY, 'genotype_data.vcf.gz'),
    GENCODE_GTF=os.path.join(os.path.dirname(__file__), EXAC_FILES_DIRECTORY, 'gencode.gtf.gz'),
    CANONICAL_TRANSCRIPT_FILE=os.path.join(os.path.dirname(__file__), EXAC_FILES_DIRECTORY, 'canonical_transcripts.txt.gz'),
    BASE_COVERAGE_FILES=[
        os.path.join(os.path.dirname(__file__), EXAC_FILES_DIRECTORY, 'coverage.txt.gz'),
    ],
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
    db.variants.ensure_index('xpos')
    db.variants.ensure_index('rsid')
    db.variants.ensure_index('genes')
    db.variants.ensure_index('transcripts')

    db.genes.remove()
    db.genes.ensure_index('gene_id')
    db.genes.ensure_index('gene_name')

    db.transcripts.remove()
    db.transcripts.ensure_index('transcript_id')
    db.transcripts.ensure_index('gene_id')

    db.exons.remove()
    db.exons.ensure_index('exon_id')
    db.exons.ensure_index('transcript_id')
    db.exons.ensure_index('gene_id')

    db.base_coverage.remove()
    db.base_coverage.ensure_index('xpos')

    # load coverage first; variant info will depend on coverage
    for filepath in app.config['BASE_COVERAGE_FILES']:
        size = os.path.getsize(filepath)
        progress = xbrowse.utils.get_progressbar(size, 'Parsing coverage: {}'.format(os.path.basename(filepath)))
        coverage_file = gzip.open(filepath)
        for base_coverage in get_base_coverage_from_file(coverage_file):
            progress.update(coverage_file.fileobj.tell())
            db.base_coverage.insert(base_coverage)
        progress.finish()

    # grab variants from sites VCF
    sites_vcf = gzip.open(app.config['SITES_VCF'])
    size = os.path.getsize(app.config['SITES_VCF'])
    progress = xbrowse.utils.get_progressbar(size, 'Loading Variants')
    for variant in get_variants_from_sites_vcf(sites_vcf):
        db.variants.insert(variant)
        progress.update(sites_vcf.fileobj.tell())
    progress.finish()

    # grab genes from GTF
    gtf_file = gzip.open(app.config['GENCODE_GTF'])
    size = os.path.getsize(app.config['GENCODE_GTF'])
    progress = xbrowse.utils.get_progressbar(size, 'Loading Genes')
    for gene in get_genes_from_gencode_gtf(gtf_file):
        db.genes.insert(gene)
        progress.update(gtf_file.fileobj.tell())
    progress.finish()
    gtf_file.close()

    # and now transcripts
    gtf_file = gzip.open(app.config['GENCODE_GTF'])
    size = os.path.getsize(app.config['GENCODE_GTF'])
    progress = xbrowse.utils.get_progressbar(size, 'Loading Transcripts')
    for transcript in get_transcripts_from_gencode_gtf(gtf_file):
        db.transcripts.insert(transcript)
        progress.update(gtf_file.fileobj.tell())
    gtf_file.close()
    progress.finish()

    # Building up gene definitions
    gtf_file = gzip.open(app.config['GENCODE_GTF'])
    size = os.path.getsize(app.config['GENCODE_GTF'])
    progress = xbrowse.utils.get_progressbar(size, 'Loading Exons')
    for exon in get_exons_from_gencode_gtf(gtf_file):
        db.exons.insert(exon)
        progress.update(gtf_file.fileobj.tell())
    gtf_file.close()
    progress.finish()

    canonical_transcript_file = gzip.open(app.config['CANONICAL_TRANSCRIPT_FILE'])
    size = os.path.getsize(app.config['CANONICAL_TRANSCRIPT_FILE'])
    progress = xbrowse.utils.get_progressbar(size, 'Loading Canonical Transcripts')
    for gene, transcript in get_canonical_transcripts(canonical_transcript_file):
        gene = db.genes.find_one({
            'gene_id': gene
        })
        if not gene:
            continue
        gene['canonical_transcript'] = transcript
        db.genes.save(gene)
        progress.update(canonical_transcript_file.fileobj.tell())
    gtf_file.close()
    progress.finish()



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
    elif datatype == 'transcript':
        return redirect('/transcript/{}'.format(identifier))
    elif datatype == 'variant':
        return redirect('/variant/{}'.format(identifier))
    elif datatype == 'region':
        return redirect('/region/{}'.format(identifier))
    elif datatype == 'dbsnp_variant_set':
        return redirect('/dbsnp/{}'.format(identifier))
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
    variant = order_variant_by_csq(variant)
    consequences = None
    if 'vep_annotations' in variant:
        consequences = defaultdict(lambda: defaultdict(list))
        for annotation in variant['vep_annotations']:
            consequences[csq_max_vep(annotation['Consequence'])][annotation['Gene']].append(annotation)
    base_coverage = lookups.get_coverage_for_bases(db, xpos, xpos + len(ref) - 1)
    any_covered = any([x['has_coverage'] for x in base_coverage])
    return render_template('variant.html', variant=variant, base_coverage=base_coverage, consequences=consequences, any_covered=any_covered)


@app.route('/gene/<gene_id>')
def gene_page(gene_id):
    db = get_db()
    gene = lookups.get_gene(db, gene_id)
    variants_in_gene = lookups.get_variants_in_gene(db, gene_id)
    transcripts_in_gene = lookups.get_transcripts_in_gene(db, gene_id)

    # Get csome anonical transcript and corresponding info
    transcript_id = gene['canonical_transcript']
    transcript = lookups.get_transcript(db, transcript_id)
    variants_in_transcript = lookups.get_variants_in_transcript(db, transcript_id)
    exons = lookups.get_exons_in_transcript(db, transcript_id)
    genomic_coord_to_exon = dict([(y, i) for i, x in enumerate(exons) for y in range(x['start'], x['stop'] + 1)])
    coverage_stats = lookups.get_coverage_for_transcript(db, genomic_coord_to_exon, transcript['xstart'], transcript['xstop'])
    add_transcript_coordinate_to_variants(db, variants_in_transcript, transcript_id)
    add_consequence_to_variants(variants_in_transcript)

    lof_variants = [
        x for x in variants_in_gene
        if any([y['LoF'] == 'HC' for y in x['vep_annotations'] if y['Gene'] == gene_id])
    ]
    composite_lof_frequency = sum([x['allele_freq'] for x in lof_variants])

    return render_template(
        'gene.html',
        gene=gene,
        transcript=transcript,
        variants_in_gene=variants_in_transcript,
        number_variants_in_gene=len(variants_in_gene),
        lof_variants_in_gene=lof_variants,
        composite_lof_frequency=composite_lof_frequency,
        transcripts_in_gene=transcripts_in_gene,
        coverage_stats=coverage_stats,
        exons=exons
    )


@app.route('/transcript/<transcript_id>')
def transcript_page(transcript_id):
    db = get_db()
    transcript = lookups.get_transcript(db, transcript_id)
    gene = lookups.get_gene(db, transcript['gene_id'])
    variants_in_transcript = lookups.get_variants_in_transcript(db, transcript_id)
    exons = lookups.get_exons_in_transcript(db, transcript_id)
    genomic_coord_to_exon = dict([(y, i) for i, x in enumerate(exons) for y in range(x['start'], x['stop'] + 1)])

    coverage_stats = lookups.get_coverage_for_transcript(db, genomic_coord_to_exon, transcript['xstart'], transcript['xstop'])

    lof_variants = [
        x for x in variants_in_transcript
        if any([y['LoF'] == 'HC' for y in x['vep_annotations'] if y['Feature'] == transcript_id])
    ]
    composite_lof_frequency = sum([x['allele_freq'] for x in lof_variants])

    add_transcript_coordinate_to_variants(db, variants_in_transcript, transcript_id)
    add_consequence_to_variants(variants_in_transcript)

    return render_template(
        'transcript.html',
        transcript=transcript,
        variants_in_transcript=variants_in_transcript,
        lof_variants=lof_variants,
        composite_lof_frequency=composite_lof_frequency,
        coverage_stats=coverage_stats,
        exons=exons,
        gene=gene
    )


@app.route('/region/<region_id>')
def region_page(region_id):
    db = get_db()
    region = region_id.split('-')
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
    return render_template(
        'region.html',
        genes_in_region=genes_in_region,
        variants_in_region=variants_in_region,
        chrom=chrom,
        start=start,
        stop=stop,
        coverage=coverage_array
    )


@app.route('/dbsnp/<rsid>')
def dbsnp_page(rsid):
    db = get_db()
    variant = lookups.get_variants_by_rsid(db, rsid)
    print variant
    return render_template(
        'region.html',
        rsid=rsid,
        variants_in_region=variant,
        coverage=None,
        genes_in_region=None
    )

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
