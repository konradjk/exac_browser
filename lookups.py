import re
from xbrowse import get_xpos


def get_gene(db, gene_id):
    return db.genes.find_one({'gene_id': gene_id}, fields={'_id': False})


def get_gene_by_name(db, gene_name):
    return db.genes.find_one({'gene_name': gene_name}, fields={'_id': False})


def get_transcript(db, transcript_id):
    return db.transcripts.find_one({'transcript_id': transcript_id}, fields={'_id': False})


def get_variant(db, xpos, ref, alt):
    return db.variants.find_one({'xpos': xpos, 'ref': ref, 'alt': alt}, fields={'_id': False})


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
            ret.append({'xpos': i})
    for item in ret:
        item['has_coverage'] = 'mean' in item
    return ret


def get_awesomebar_suggestions(db, query):
    """
    This generates autocomplete suggestions when user
    query is the string that user types
    If it is the prefix for a gene, return list of gene names
    """
    regex = re.compile('^' + re.escape(query), re.IGNORECASE)
    print regex
    genes = db.genes.find({'gene_name': {
        '$regex': regex,
    }}).limit(20)
    genes = list(genes)
    if genes is None:
        genes = []
    return [gene['gene_name'] for gene in genes]


# 1:1-1000
R1 = re.compile(r'^(\d+|X|Y|M|MT):(\d+)-(\d+)$')

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
    variant = db.variants.find_one({'rsid' : query}) # TODO - https://github.com/brettpthomas/exac_browser/issues/19
    if variant:
        # TODO - https://github.com/brettpthomas/exac_browser/issues/18
        return 'variant', variant
    # TODO - https://github.com/brettpthomas/exac_browser/issues/14

    # Region
    m = R1.match(query)
    if m:
        return 'region', '{}-{}-{}'.format(m.group(1), m.group(2), m.group(3))

    


def get_genes_in_region(db, chrom, start, stop):
    """
    Genes that overlap a region
    """
    xstart = get_xpos(chrom, start)
    xstop = get_xpos(chrom, stop)
    genes = db.genes.find({
        'xstart': {'$lte': xstop},
        'xstop': {'$gte': xstart},
    })
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
    })
    return list(variants)


def get_variants_in_gene(db, gene_id):
    """
    """
    return list(db.variants.find({'genes': gene_id}, fields={'_id': False}))


def get_transcripts_in_gene(db, gene_id):
    """
    """
    return list(db.transcripts.find({'gene_id': gene_id}, fields={'_id': False}))


def get_variants_in_transcript(db, transcript_id):
    """
    """
    return list(db.variants.find({'transcripts': transcript_id}, fields={'_id': False}))
