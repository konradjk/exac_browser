import re
import itertools
from xbrowse import get_xpos
from utils import *

SEARCH_LIMIT = 10000
UNSUPPORTED_QUERIES = ['TTN', 'ENSG00000155657', 'CMD1G', 'CMH9', 'CMPD4', 'FLJ32040', 'LGMD2J', 'MYLK5', 'TMD',
                       u'ENST00000342175', u'ENST00000359218', u'ENST00000342992', u'ENST00000460472',
                       u'ENST00000589042', u'ENST00000591111']


def get_gene(db, gene_id):
    return db.genes.find_one({'gene_id': gene_id}, fields={'_id': False})


def get_gene_by_name(db, gene_name):
    # try gene_name field first
    gene = db.genes.find_one({'gene_name': gene_name}, fields={'_id': False})
    if gene:
        return gene
    # if not, try gene['other_names']
    return db.genes.find_one({'other_names': gene_name}, fields={'_id': False})


def get_transcript(db, transcript_id):
    transcript = db.transcripts.find_one({'transcript_id': transcript_id}, fields={'_id': False})
    if not transcript:
        return None
    transcript['exons'] = get_exons_in_transcript(db, transcript_id)
    return transcript


def get_variant(db, xpos, ref, alt):
    return db.variants.find_one({'xpos': xpos, 'ref': ref, 'alt': alt}, fields={'_id': False})


def get_variants_by_rsid(db, rsid):
    if not rsid.startswith('rs'):
        return None
    try:
        int(rsid.lstrip('rs'))
    except Exception, e:
        return None
    variants = list(db.variants.find({'rsid': rsid}, fields={'_id': False}))
    add_consequence_to_variants(variants)
    return variants


def get_coverage_for_bases(db, xstart, xstop=None):
    """
    Get the coverage for the list of bases given by xstart->xstop, inclusive
    Returns list of coverage dicts
    xstop can be None if just one base, but you'll still get back a list
    """
    if xstop is None:
        xstop = xstart
    coverages = {
        doc['xpos']: doc for doc in db.base_coverage.find(
            {'xpos': {'$gte': xstart, '$lte': xstop}},
            fields={'_id': False}
        )
    }
    ret = []
    for i in range(xstart, xstop+1):
        if i in coverages:
            ret.append(coverages[i])
        else:
            ret.append({'xpos': i, 'pos': xpos_to_pos(i)})
    for item in ret:
        item['has_coverage'] = 'mean' in item
    return ret


def get_coverage_for_transcript(db, xstart, xstop=None):
    """

    :param db:
    :param genomic_coord_to_exon:
    :param xstart:
    :param xstop:
    :return:
    """
    coverage_array = get_coverage_for_bases(db, xstart, xstop)
    # only return coverages that have coverage (if that makes any sense?)
    return coverage_array
    #return [c for c in coverage_array if c['has_coverage']]


def get_awesomebar_suggestions(db, query):
    """
    This generates autocomplete suggestions when user
    query is the string that user types
    If it is the prefix for a gene, return list of gene names
    """
    regex = re.compile('^' + re.escape(query), re.IGNORECASE)
    genes = db.genes.find({'gene_name': {
        '$regex': regex,
    }}, fields={'gene_name': True}).limit(20)
    if genes is None:
        genes = []
    gene_names = sorted([gene['gene_name'] for gene in genes])

    genes_by_other_name = db.genes.find({'other_names': {
        '$regex': regex,
    }}, fields={'other_names': True}).limit(20)
    if genes_by_other_name is None:
        genes_by_other_name = []
    other_names = [gene['other_names'] for gene in genes_by_other_name]
    other_names = sorted([name for names in other_names for name in names if regex.match(name)])
    return gene_names + other_names


# 1:1-1000
R1 = re.compile(r'^(\d+|X|Y|M|MT)\s*:\s*(\d+)-(\d+)$')
R2 = re.compile(r'^(\d+|X|Y|M|MT)\s*:\s*(\d+)$')
R3 = re.compile(r'^(\d+|X|Y|M|MT)$')
R4 = re.compile(r'^(\d+|X|Y|M|MT)\s*[-:]\s*(\d+)-([ATCG]+)-([ATCG]+)$')


def get_awesomebar_result(db, query):
    """
    Similar to the above, but this is after a user types enter
    We need to figure out what they meant - could be gene, variant, region

    Return tuple of (datatype, identifier)
    Where datatype is one of 'gene', 'variant', or 'region'
    And identifier is one of:
    - ensembl ID for gene
    - variant ID string for variant (eg. 1-1000-A-T)
    - region ID string for region (eg. 1-1000-2000)

    Follow these steps:
    - if query is an ensembl ID, return it
    - if a gene symbol, return that gene's ensembl ID
    - if an RSID, return that variant's string


    Finally, note that we don't return the whole object here - only it's identifier.
    This could be important for performance later

    """
    query = query.strip()
    print query
    if query in UNSUPPORTED_QUERIES:
        return 'error', query

    # Variant
    variant = get_variants_by_rsid(db, query.lower())
    if variant:
        if len(variant) == 1:
            return 'variant', variant[0]['variant_id']
        else:
            return 'dbsnp_variant_set', variant[0]['rsid']
    # variant = get_variant(db, )
    # TODO - https://github.com/brettpthomas/exac_browser/issues/14

    # From here out, all should be uppercase (gene, tx, region, variant_id)
    query = query.upper()

    # Ensembl formatted queries
    if query.startswith('ENS'):
        # Gene
        gene = get_gene(db, query)
        if gene:
            return 'gene', gene['gene_id']

        # Transcript
        transcript = get_transcript(db, query)
        if transcript:
            return 'transcript', transcript['transcript_id']

    gene = get_gene_by_name(db, query)
    if gene:
        return 'gene', gene['gene_id']

    # From here on out, only region queries
    query = query.lstrip('chr')
    # Region
    m = R1.match(query)
    if m:
        if int(m.group(3)) < int(m.group(2)):
            return 'region', 'invalid'
        return 'region', '{}-{}-{}'.format(m.group(1), m.group(2), m.group(3))
    m = R2.match(query)
    if m:
        return 'region', '{}-{}-{}'.format(m.group(1), m.group(2), m.group(2))
    m = R3.match(query)
    if m:
        return 'region', '{}'.format(m.group(1))
    m = R4.match(query)
    if m:
        return 'variant', '{}-{}-{}-{}'.format(m.group(1), m.group(2), m.group(3), m.group(4))

    return 'error', query


def get_genes_in_region(db, chrom, start, stop):
    """
    Genes that overlap a region
    """
    xstart = get_xpos(chrom, start)
    xstop = get_xpos(chrom, stop)
    genes = db.genes.find({
        'xstart': {'$lte': xstop},
        'xstop': {'$gte': xstart},
    }, fields={'_id': False})
    return list(genes)


def get_variants_in_region(db, chrom, start, stop):
    """
    Variants that overlap a region
    Unclear if this will include CNVs
    """
    xstart = get_xpos(chrom, start)
    xstop = get_xpos(chrom, stop)
    variants = db.variants.find({
        'xstart': {'$lte': xstop},  # start of variant should be before (or equal to) end of region
        'xstop': {'$gte': xstart},  # opposite of above
    }, fields={'_id': False}, limit=SEARCH_LIMIT)
    add_consequence_to_variants(variants)
    return list(variants)


def remove_extraneous_information(variant):
    del variant['genotype_depths']
    del variant['genotype_qualities']
    del variant['pop_acs']
    del variant['pop_ans']
    variant['hom_count'] = sum(variant['pop_homs'].values())
    del variant['pop_homs']
    del variant['transcripts']
    del variant['genes']
    del variant['orig_alt_alleles']
    del variant['xpos']
    del variant['xstart']
    del variant['xstop']
    del variant['site_quality']
    del variant['vep_annotations']


def get_variants_in_gene(db, gene_id):
    """
    """
    variants = []
    for variant in db.variants.find({'genes': gene_id}, fields={'_id': False}):
        add_consequence_to_variant(variant)
        remove_extraneous_information(variant)
        variants.append(variant)
    return variants


def get_transcripts_in_gene(db, gene_id):
    """
    """
    return list(db.transcripts.find({'gene_id': gene_id}, fields={'_id': False}))


def get_variants_in_transcript(db, transcript_id):
    """
    """
    variants = []
    for variant in db.variants.find({'transcripts': transcript_id}, fields={'_id': False}):
        variant['vep_annotations'] = [x for x in variant['vep_annotations'] if x['Feature'] == transcript_id]
        add_consequence_to_variant(variant)
        remove_extraneous_information(variant)
        variants.append(variant)
    return variants


def get_exons_in_transcript(db, transcript_id):
    # return sorted(
    #     [x for x in
    #      db.exons.find({'transcript_id': transcript_id}, fields={'_id': False})
    #      if x['feature_type'] != 'exon'],
    #     key=lambda k: k['start'])
    return sorted(list(db.exons.find({'transcript_id': transcript_id, 'feature_type': { "$in": ['CDS', 'UTR', 'exon'] }}, fields={'_id': False})), key=lambda k: k['start'])
