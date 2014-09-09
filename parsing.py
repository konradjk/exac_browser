"""
Utils for reading flat files that are loaded into database
"""
import re
import copy
import xbrowse


def get_base_coverage_from_file(base_coverage_file):
    """
    Read a base coverage file and return iter of dicts that look like:
    {
        'xpos': 1e9+1,
        'mean': 0.0,
        'median': 0.0,
        '1': 0.0,
        '5': 0.0,
        '10': 0.0,
        '15': 0.0,
        '20': 0.0,
        '25': 0.0,
        '30': 0.0,
        '50': 0.0,
        '100': 0.0,
    }
    """
    pass


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
        info_field = dict([(x.split('=', 1)) if '=' in x else (x, x) for x in re.split(';(?=\w)', fields[7])])
        consequence_array = info_field['CSQ'].split(',') if 'CSQ' in info_field else []
        annotations = [dict(zip(vep_field_names, x.split('|'))) for x in consequence_array if len(vep_field_names) == len(x.split('|'))]

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
            variant['xstart'] = variant['xpos']
            variant['xstop'] = variant['xpos'] + len(variant['alt']) - len(variant['ref'])
            variant['variant_id'] = '{}-{}-{}-{}'.format(variant['chrom'], str(variant['pos']), variant['ref'], variant['alt'])
            variant['orig_alt_alleles'] = alt_alleles
            variant['site_quality'] = float(fields[5])
            variant['filter'] = fields[6]
            variant['vep_annotations'] = [ann for ann in annotations if int(ann['ALLELE_NUM']) == i + 1]
            variant['allele_count'] = int(info_field['AC'].split(',')[i])
            variant['allele_freq'] = float(info_field['AF'].split(',')[i])
            variant['num_alleles'] = int(info_field['AN'])
            variant['genes'] = list({annotation['Gene'] for annotation in variant['vep_annotations']})
            variant['transcripts'] = list({annotation['Feature'] for annotation in variant['vep_annotations']})

            yield variant


def get_genotype_data_from_full_vcf(full_vcf):
    """

    """
    for line in full_vcf:
        line = line.strip('\n')
        if line.startswith('#'):
            continue
        fields = line.split('\t')
        info_field = dict([(x.split('=', 1)) if '=' in x else (x, x) for x in re.split(';(?=\w)', fields[7])])

        alt_alleles = fields[4].split(',')

        format_field = fields[8].split(':')
        format_data = [dict(zip(format_field, x.split(':'))) for x in fields[9:]]

        # different variant for each alt allele
        for i, alt_allele in enumerate(alt_alleles):

            genotype_info_container = {
                'xpos': xbrowse.get_xpos(fields[0], int(fields[1])),
                'ref': fields[3],
                'alt': alt_allele,
                'genotype_info': {
                    'something': 0.4,
                    'genotype_qualities': [int(x['GQ']) for x in format_data if 'GQ' in x],
                    # I don't love this: I thought DP had to be int (even if 0), but apparently we have '.'
                    'genotype_depths': [int(x['DP']) if x['DP'].isdigit() else 0 for x in format_data if 'DP' in x],
                    'genotypes': [re.split("/|\|", x['GT']) for x in format_data if 'GT' in x],
                }
            }
            yield genotype_info_container

def get_genes_from_gencode_gtf(gtf_file):
    """
    Parse gencode GTF file;
    Returns iter of gene dicts
    """
    for line in gtf_file:
        if line.startswith('#'):
            continue
        fields = line.strip('\n').split('\t')

        # This isn't required anymore, as we're currently reading in just a genes GTF, but could be useful later
        if fields[2] != 'gene':
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


def get_transcripts_from_gencode_gtf(gtf_file):
    """
    Parse gencode GTF file;
    Returns iter of transcript dicts
    """
    for line in gtf_file:
        if line.startswith('#'):
            continue
        fields = line.strip('\n').split('\t')

        if fields[2] != 'transcript':
            continue

        chrom = fields[0][3:]
        start = int(fields[3]) + 1  # bed files are 0-indexed
        stop = int(fields[4]) + 1
        info = dict(x.strip().split() for x in fields[8].split(';') if x != '')
        info = {k: v.strip('"') for k, v in info.items()}
        transcript_id = info['transcript_id'].split('.')[0]
        gene_id = info['gene_id'].split('.')[0]

        gene = {
            'transcript_id': transcript_id,
            'gene_id': gene_id,
            'chrom': chrom,
            'start': start,
            'stop': stop,
            'xstart': xbrowse.get_xpos(chrom, start),
            'xstop': xbrowse.get_xpos(chrom, stop),
        }
        yield gene