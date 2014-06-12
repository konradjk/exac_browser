

def get_gene(db, gene_id):
    return db.genes.find_one({'gene_id': gene_id})


def get_variant(db, xpos, ref, alt):
    return db.variants.find_one({'xpos': xpos, 'ref': ref, 'alt': alt})


def get_awesomebar_suggestions(db, query):
    """
    This generates autocomplete suggestions when user
    query is the string that user types
    If it is the prefix for a gene, return list of gene names
    QUESTION: should we autocomplete anything else?
    """
    raise NotImplementedError


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
    1. if query is an ensembl ID, return it
    2. if a gene symbol, return that gene's ensembl ID
    3. if an RSID, return that variant's string

    Finally, note that we don't return the whole object here - only it's identifier.
    This could be important for performance later

    """
    gene = db.genes.find_one({'gene_id': query})
    if gene:
        return 'gene', gene['gene_id']
    gene = db.genes.find_one({'gene_name': query})
    if gene:
        return 'gene', gene['gene_id']



def get_genes_in_region(db, xstart, xstop):
    """
    Genes that overlap a region
    """
    raise NotImplementedError


def get_variants_in_region(db, xstart, xstop):
    """
    Variants that overlap a region
    Unclear if this will include CNVs
    """
    raise NotImplementedError