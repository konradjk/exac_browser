"""
Utils for reading flat files that are loaded into database
"""
import re
import copy
import xbrowse


def get_variants_from_sites_vcf(sites_vcf):
    """
    Parse exac sites VCF file and return iter of variant dicts
    sites_vcf is a file (gzipped), not file path
    """
    vep_field_names = None
    for line in sites_vcf:
        line = line.strip('\n')
        if line.startswith('##INFO=<ID=CSQ'):
            vep_field_names = line.split('Format: ')[-1].strip('">').split('|')
        if line.startswith('#'):
            continue

        # If we get here, it's a variant line

        # This elegant parsing code below is copied from https://github.com/konradjk/loftee
        fields = line.split('\t')
        info_field = dict([(x.split('=', 1)) for x in re.split(';(?=\w)', fields[7]) if x.find('=') > -1])
        if 'CSQ' not in info_field:
            info_field['CSQ'] = ''
        annotations = [dict(zip(vep_field_names, x.split('|'))) for x in info_field['CSQ'].split(',')]

        alt_alleles = fields[4].split(',')

        # different variant for each alt allele
        for i, alt_allele in enumerate(alt_alleles):

            # Variant is just a dict
            # Make a copy of the info_field dict - so all the original data remains
            # Add some new keys that are allele-specific
            variant = {}
            variant['chrom'] = fields[0]
            variant['pos'] = int(fields[1])
            variant['rsid'] = fields[2]
            variant['xpos'] = xbrowse.get_xpos(variant['chrom'], variant['pos'])
            variant['ref'] = fields[3]
            variant['alt'] = alt_allele
            variant['orig_alt_alleles'] = alt_alleles
            variant['site_quality'] = float(fields[5])
            variant['filter'] = fields[6]
            variant['vep_annotations'] = [ann for ann in annotations if int(ann['ALLELE_NUM']) == i]
            variant['allele_count'] = int(info_field['AC'].split(',')[i])
            variant['allele_freq'] = float(info_field['AF'].split(',')[i])
            variant['num_alleles'] = int(info_field['AN'])

            yield variant


def get_genes_from_gencode_gtf(gtf_file):
    """
    Parse gencode GTF file;
    Returns iter of gene dicts
    """
    for line in gtf_file:
        if line.startswith('#'):
            continue
        fields = line.strip('\n').split('\t')

        # only look at ensembl genes. may want to change this
        if fields[1] != 'ENSEMBL' and fields[2] != 'gene':
            continue

        chrom = fields[0][3:]
        start = int(fields[3]) + 1  # bed files are 0-indexed
        stop = int(fields[4]) + 1
        info = dict(x.strip().split() for x in fields[8].split(';') if x != '')
        info = {k: v.strip('"') for k, v in info.items()}
        gene_id = info['gene_id'].split('.')[0]

        gene = {
            'gene_id': gene_id,
            'gene_name': info['gene_name'],
            'chrom': chrom,
            'start': start,
            'stop': stop,
            'xstart': xbrowse.get_xpos(chrom, start),
            'xstop': xbrowse.get_xpos(chrom, stop),
        }
        yield gene
