from operator import itemgetter


def add_transcript_coordinate_to_variants(db, variant_list, transcript_id):
    """
    Each variant has a 'xpos' and 'pos' positional attributes.
    This method takes a list of variants and adds a third position: the "transcript coordinates".
    This is defined as the distance from the start of the transcript, in coding bases.
    So a variant in the 7th base of the 6th exon of a transcript will have a transcript coordinate of
    the sum of the size of the first 5 exons) + 7
    This is 0-based, so a variant in the first base of the first exon has a transcript coordinate of 0.

    You may want to add transcript coordinates for multiple transcripts, so this is stored in a variant as
    variant['transcript_coordinates'][transcript_id]

    If a variant in variant_list does not have a `transcript_coordinates` dictionary, we create one

    If a variant start position for some reason does not fall in any exons in this transcript, its coordinate is 0.
    This is perhaps logically inconsistent,
    but it allows you to spot errors quickly if there's a pileup at the first base.
    `None` would just break things.

    Consider the behavior if a 20 base deletion deletes parts of two exons.
    I think the behavior in this method is consistent, but beware that it might break things downstream.

    Edits variant_list in place; no return val
    """

    import lookups
    # make sure exons is sorted by (start, end)
    exons = sorted(lookups.get_exons_in_transcript(db, transcript_id), key=itemgetter('start', 'stop'))

    # offset from start of base for exon in ith position (so first item in this list is always 0)
    exon_offsets = [0 for i in range(len(exons))]
    for i, exon in enumerate(exons):
        for j in range(i+1, len(exons)):
            exon_offsets[j] += exon['stop'] - exon['start']

    for variant in variant_list:
        if 'transcript_coordinates' not in variant:
            variant['transcript_coordinates'] = {}
        variant['transcript_coordinates'][transcript_id] = 0
        for i, exon in enumerate(exons):
            if exon['start'] <= variant['pos'] <= exon['stop']:
                variant['transcript_coordinates'][transcript_id] = exon_offsets[i] + variant['pos'] - exon['start']


def xpos_to_pos(xpos):
    return int(xpos % 1e9)


def add_consequence_to_variants(variant_list):
    for variant in variant_list:
        variant['major_consequence'] = csq_max([csq_max_vep(x['Consequence']) for x in variant['vep_annotations']])
        if csq_order_dict[variant['major_consequence']] <= csq_order_dict["initiator_codon_variant"]:
            variant['category'] = 'lof_variant'
        elif csq_order_dict[variant['major_consequence']] <= csq_order_dict["missense_variant"]:
            # Should be noted that this grabs inframe deletion, etc.
            variant['category'] = 'missense_variant'
        else:
            variant['category'] = 'other_variant'


csq_order = ["transcript_ablation",
"splice_donor_variant",
"splice_acceptor_variant",
"stop_gained",
"frameshift_variant",
"stop_lost",
"initiator_codon_variant",
"inframe_insertion",
"inframe_deletion",
"missense_variant",
"transcript_amplification",
"splice_region_variant",
"incomplete_terminal_codon_variant",
"synonymous_variant",
"stop_retained_variant",
"coding_sequence_variant",
"mature_miRNA_variant",
"5_prime_UTR_variant",
"3_prime_UTR_variant",
"non_coding_exon_variant",
"nc_transcript_variant",
"intron_variant",
"NMD_transcript_variant",
"upstream_gene_variant",
"downstream_gene_variant",
"TFBS_ablation",
"TFBS_amplification",
"TF_binding_site_variant",
"regulatory_region_variant",
"regulatory_region_ablation",
"regulatory_region_amplification",
"feature_elongation",
"feature_truncation",
"intergenic_variant",
""]
csq_order_dict = dict(zip(csq_order, range(len(csq_order))))
rev_csq_order_dict = dict(zip(range(len(csq_order)), csq_order))


def csq_max_vep(ann_list):
    return rev_csq_order_dict[csq_max_list(ann_list.split('&'))]


def csq_max(ann_list):
    if len(ann_list) == 0: return ''
    return rev_csq_order_dict[csq_max_list(ann_list)]


def csq_max_list(ann_list):
    return min([csq_order_dict[ann] for ann in ann_list])


def order_variant_by_csq(annotation_list):
    # new_order = [csq_order_dict[ann] for ann in annotation_list]
    # TODO: get order
    return annotation_list
    #return [annotation_list[x] for x in new_order]


def get_minimal_representation(pos, ref, alt): 
    """
    Get the minimal representation of a variant, based on the ref + alt alleles in a VCF
    This is used to make sure that multiallelic variants in different datasets, 
    with different combinations of alternate alleles, can always be matched directly. 

    Note that chromosome is ignored here - in xbrowse, we'll probably be dealing with 1D coordinates 
    Args: 
        pos (int): genomic position in a chromosome (1-based)
        ref (str): ref allele string
        alt (str): alt allele string
    Returns: 
        tuple: (pos, ref, alt) of remapped coordinate
    """
    pos = int(pos)
    # If it's a simple SNV, don't remap anything
    if len(ref) == 1 and len(alt) == 1: 
        return pos, ref, alt
    else:
        # strip off identical suffixes
        while(alt[-1] == ref[-1] and min(len(alt),len(ref)) > 1):
            alt = alt[:-1]
            ref = ref[:-1]
        # strip off identical prefixes and increment position
        while(alt[0] == ref[0] and min(len(alt),len(ref)) > 1):
            alt = alt[1:]
            ref = ref[1:]
            pos += 1
        return pos, ref, alt 
