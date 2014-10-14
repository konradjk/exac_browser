import re
from xbrowse import get_xpos
from utils import xpos_to_pos

SEARCH_LIMIT = 10000

def get_gene(db, gene_id):
    return db.genes.find_one({'gene_id': gene_id}, fields={'_id': False})


def get_gene_by_name(db, gene_name):
    return db.genes.find_one({'gene_name': gene_name}, fields={'_id': False})


def get_transcript(db, transcript_id):
    transcript = db.transcripts.find_one({'transcript_id': transcript_id}, fields={'_id': False})
    if not transcript:
        return None
    transcript['exons'] = get_exons_in_transcript(db, transcript_id)
    return transcript


def get_variant(db, xpos, ref, alt):
    return db.variants.find_one({'xpos': xpos, 'ref': ref, 'alt': alt}, fields={'_id': False})


def get_variants_by_rsid(db, rsid):
    print 'Looking up: ', rsid
    if not rsid.startswith('RS'):
        return None
    return list(db.variants.find({'rsid': rsid}, fields={'_id': False}))


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
    }}).limit(20)
    genes = list(genes)
    if genes is None:
        genes = []
    return [gene['gene_name'] for gene in genes]


# 1:1-1000
R1 = re.compile(r'^(\d+|X|Y|M|MT)\s*:\s*(\d+)-(\d+)$')
R2 = re.compile(r'^(\d+|X|Y|M|MT)\s*:\s*(\d+)$')
R3 = re.compile(r'^(\d+|X|Y|M|MT)$')
# R1 = re.compile(r'^(\d+|X|Y|M|MT):(\d+)-(\d+)$')
# R2 = re.compile(r'^(\d+|X|Y|M|MT):(\d+)$')

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
    query = query.strip().upper()
    print query

    # Gene
    gene = get_gene(db, query)
    if gene:
        return 'gene', gene['gene_id']
    gene = get_gene_by_name(db, query)
    if gene:
        return 'gene', gene['gene_id']

    # Transcript
    transcript = get_transcript(db, query)
    if transcript:
        return 'transcript', transcript['transcript_id']

    # Variant
    variant = get_variants_by_rsid(db, query)
    if variant:
        if len(variant) == 1:
            return 'variant', variant[0]['variant_id']
        else:
            print 'Got ', variant
            return 'dbsnp_variant_set', variant[0]['rsid']
    # variant = get_variant(db, )
    # TODO - https://github.com/brettpthomas/exac_browser/issues/14

    # Region
    m = R1.match(query.lstrip('chr'))
    if m:
        if int(m.group(3)) < int(m.group(2)):
            return 'region', 'invalid'
        return 'region', '{}-{}-{}'.format(m.group(1), m.group(2), m.group(3))
    m = R2.match(query.lstrip('chr'))
    if m:
        return 'region', '{}-{}-{}'.format(m.group(1), m.group(2), m.group(2))
    m = R3.match(query.lstrip('chr'))
    if m:
        return 'region', '{}'.format(m.group(1))
    print "Didn't find anything"

    


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
    return list(variants)


def get_variants_in_gene(db, gene_id):
    """
    """
    variants = list(db.variants.find({'genes': gene_id}, fields={'_id': False}))
    for variant in variants:
        del variant['genotype_depths']
        del variant['genotype_qualities']
        del variant['pop_acs']
        del variant['pop_ans']
        del variant['pop_homs']

        csq = variant['vep_annotations']
        variant['vep_annotations'] = [
            {'Consequence': x['Consequence'],
             'Gene': x['Gene'],
             'Feature': x['Feature'],
             'LoF': x['LoF']}
            for x in csq
        ]
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
        variants.append(variant)
    return variants


def get_exons_in_transcript(db, transcript_id):
    # return sorted(
    #     [x for x in
    #      db.exons.find({'transcript_id': transcript_id}, fields={'_id': False})
    #      if x['feature_type'] != 'exon'],
    #     key=lambda k: k['start'])
    return sorted(list(db.exons.find({'transcript_id': transcript_id, 'feature_type': { "$in": ['CDS', 'UTR'] }}, fields={'_id': False})), key=lambda k: k['start'])
