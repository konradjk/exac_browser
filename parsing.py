"""
Utils for reading flat files that are loaded into database
"""
import copy
import csv
import re
import traceback

from utils import *


POPS = {
    'AFR': 'African',
    'AMR': 'Latino',
    'EAS': 'East Asian',
    'EUR': 'European',
    'FIN': 'European (Finnish)',
    'NFE': 'European (Non-Finnish)',
    'SAS': 'South Asian',
    'OTH': 'Other'
}
POP_NUM = {
    'African American': 2530,
    'European American': 7239,
}
POP_NUM_RECQL = {
    'African American': 1260,
    'European American': 3578,
}
ALLELE_NUM = sum(POP_NUM.values())
ALLELE_NUM_RECQL = sum(POP_NUM_RECQL.values())


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

    float_header_fields = ['mean', 'median', '1', '5', '10', '15', '20', '25', '30', '50', '100']
    for line in base_coverage_file:
        if line.startswith('#'):
            continue
        fields = line.strip('\n').split('\t')
        d = {
            'xpos': get_xpos(fields[0], int(fields[1])),
            'pos': int(fields[1]),
        }
        for i, k in enumerate(float_header_fields):
            d[k] = float(fields[i+2])
        yield d


def get_variants_from_sites_vcf(sites_vcf, all_transcripts, gene_ids_by_name=None):
    """
    Parse exac sites VCF file and return iter of variant dicts
    sites_vcf is a file (gzipped), not file path
    all_transcripts is a set of string transcript ids
    gene_ids_by_name is a dict mapping gene_name_upper to gene_id
    """
    if not gene_ids_by_name:
        gene_ids_by_name = {}

    vep_field_names = []
    eff_field_names = []
    lof_field_names = []
    nmd_field_names = []

    for line in sites_vcf:
        try:
            line = line.strip('\n')
            if line.startswith('##INFO=<ID=CSQ'):
                vep_field_names = line.split('Format: ')[-1].strip('">').split('|')
            if line.startswith('##INFO=<ID=EFF'):
                eff_field_names = [field.strip(' []()') for field in
                    line.split("Format: 'Effect (")[-1].strip(')\' ">').split('|')]
            if line.startswith('##INFO=<ID=LOF'):
                lof_field_names = [field.strip(' []()') for field in
                    line.split("Format: '")[-1].strip(" '>").split('|')]
            if line.startswith('##INFO=<ID=NMD'):
                nmd_field_names = [field.strip(' []()') for field in
                    line.split("Format: '")[-1].strip(" '>").split('|')]
            if line.startswith('##INFO=<ID=DP_HIST'):
                dp_mids = map(float, line.split('Mids: ')[-1].strip('">').split('|'))
            if line.startswith('##INFO=<ID=GQ_HIST'):
                gq_mids = map(float, line.split('Mids: ')[-1].strip('">').split('|'))
            if line.startswith('#'):
                continue

            # If we get here, it's a variant line
            # This elegant parsing code below is copied from https://github.com/konradjk/loftee
            fields = line.split('\t')
            info_field = dict([(x.split('=', 1)) if '=' in x else (x, x) for x in re.split(';(?=\w)', fields[7])])
            annotations = [dict(zip(vep_field_names, x.split('|'))) for x in info_field.get('CSQ', '').split(',')]
            coding_annotations = [ann for ann in annotations if ann.get('Feature', '').startswith('ENST')]

            eff_annotations = []
            transcripts = set()
            for eff in info_field.get('EFF', '').split(','):
                eff = eff.strip(' )')
                if not eff:
                    continue
                effect, rest = eff.split('(')
                values = dict(zip(eff_field_names, rest.split('|')))
                values['Effect'] = effect
                transcript = values.get('Transcript_ID')
                if transcript:
                    if transcript not in all_transcripts:
                        # this transcript isn't one of ours, so ignore the annotation
                        continue
                    transcripts.add(transcript)
                eff_annotations.append(values)

            lof_annotations = [lof.strip(' )(') for lof in info_field.get('LOF', '').split('|')]
            nmd_annotations = [nmd.strip(' )(') for nmd in info_field.get('NMD', '').split('|')]

            alt_alleles = fields[4].split(',')

            # different variant for each alt allele
            for i, alt_allele in enumerate(alt_alleles):
                vep_annotations = [ann for ann in coding_annotations if int(ann['ALLELE_NUM']) == i + 1]

                transcripts |= {a['Feature'] for a in vep_annotations}
                if not transcripts & all_transcripts:
                    # this variant isn't on any of our transcripts, so ignore it
                    continue

                # Variant is just a dict
                # Make a copy of the info_field dict - so all the original data remains
                # Add some new keys that are allele-specific
                pos, ref, alt = get_minimal_representation(fields[1], fields[3], alt_allele)

                variant = {}
                variant['chrom'] = fields[0]
                variant['pos'] = pos
                variant['rsid'] = fields[2]
                variant['xpos'] = get_xpos(variant['chrom'], variant['pos'])
                variant['ref'] = ref
                variant['alt'] = alt
                variant['xstart'] = variant['xpos']
                variant['xstop'] = variant['xpos'] + len(variant['alt']) - len(variant['ref'])
                variant['variant_id'] = '{}-{}-{}-{}'.format(variant['chrom'], variant['pos'], variant['ref'], variant['alt'])
                variant['orig_alt_alleles'] = [
                    '{}-{}-{}-{}'.format(variant['chrom'], *get_minimal_representation(fields[1], fields[3], x))
                    for x in alt_alleles
                ]
                try:
                    variant['site_quality'] = float(fields[5])
                except ValueError:
                    variant['site_quality'] = 0
                variant['filter'] = fields[6]
                if variant['filter'] == '.':
                    variant['filter'] = 'PASS'
                variant['vep_annotations'] = vep_annotations
                variant['eff_annotations'] = eff_annotations
                variant['lof_annotations'] = lof_annotations
                variant['nmd_annotations'] = nmd_annotations

                AC = info_field.get('AC_Adj') or info_field.get('AC')
                AN = info_field.get('AN_Adj') or info_field.get('AN')
                if AC and AN:
                    # TODO XXX fix clr vcf:accumulate to emit AC per alt allele,
                    # then switch this index back to i?
                    variant['allele_count'] = int(AC.split(',')[0])
                    if not variant['allele_count'] and variant['filter'] == 'PASS':
                        variant['filter'] = 'AC_Adj0' # Temporary filter
                    variant['allele_num'] = int(AN)

                    if variant['allele_num'] > 0:
                        variant['allele_freq'] = variant['allele_count'] / float(variant['allele_num'])
                    else:
                        variant['allele_freq'] = None

                for field in 'AC', 'AN', 'Hom':
                    pops = variant['pop_%ss' % field.lower()] = {}
                    for code, desc in POPS.items():
                        val = info_field.get('%s_%s' % (field, code))
                        if val:
                            pops[desc] = int(val)
                variant['hom_count'] = sum(variant.get('pop_homs', {}).values())

                if ('AC_MALE' in info_field and 'AC_FEMALE' in info_field and
                    'AN_MALE' in info_field and 'AN_FEMALE' in info_field):
                    variant['ac_male'] = info_field['AC_MALE']
                    variant['ac_female'] = info_field['AC_FEMALE']
                    variant['an_male'] = info_field['AN_MALE']
                    variant['an_female'] = info_field['AN_FEMALE']

                if variant['chrom'] in ('X', 'Y') and 'Hemi_%s' % pop in info_field:
                    variant['pop_hemis'] = {desc: int(info_field['Hemi_%s' % x].split(',')[i]) for code, desc in POPS.items()}
                    variant['hemi_count'] = sum(variant['pop_hemis'].values())

                variant['quality_metrics'] = {x: info_field[x] for x in METRICS if x in info_field}

                variant['genes'] = list(
                    {a['Gene'] for a in vep_annotations} |
                    {gene_ids_by_name[a['Gene_Name'].upper()] for a in eff_annotations if a.get('Gene_Name')})
                variant['transcripts'] = list(transcripts)

                if 'DP_HIST' in info_field:
                    hists_all = [info_field['DP_HIST'].split(',')[0], info_field['DP_HIST'].split(',')[i+1]]
                    variant['genotype_depths'] = [zip(dp_mids, map(int, x.split('|'))) for x in hists_all]
                if 'GQ_HIST' in info_field:
                    hists_all = [info_field['GQ_HIST'].split(',')[0], info_field['GQ_HIST'].split(',')[i+1]]
                    variant['genotype_qualities'] = [zip(gq_mids, map(int, x.split('|'))) for x in hists_all]

                yield variant
        except Exception:
            print("Error parsing vcf line: " + line)
            traceback.print_exc()
            raise


def get_variants_from_whi_tsv(tsv_file, all_transcripts, gene_ids_by_name=None):
    """
    Parse a TSV variant file for the WHI project, returns iter of variant dicts.

    sites_tsv is a file object, not file path
    all_transcripts is a set of string transcript ids
    gene_ids_by_name is a dict mapping gene_name_upper to gene_id
    """
    if not gene_ids_by_name:
        gene_ids_by_name = {}

    types = {  # maps line['Type'] to utils.csq_order
        'UTR3': '3_prime_UTR_variant',
        'UTR5': '5_prime_UTR_variant',
        'delFS': 'frameshift_deletion',
        'delIF': 'inframe_deletion',
        'insFS': 'frameshift_insertion',
        'insIF': 'inframe_insertion',
        'intronic': 'intron_variant',
        'missense': 'missense_variant',
        'stoploss': 'stop_lost',
        # 'nonsense': '',
        # 'silent': '',
        # 'splice': '',
    }

    for line in csv.DictReader(tsv_file, dialect='excel-tab'):
        # normalize '-' to blank or 0
        numeric = {'Start', 'All carrier count', 'Homoz', 'AA', 'EA',
                   'Max_qual', 'Max_MQ'}
        line = {k: (0 if k in numeric else '') if v == '-' else v
                for k, v in line.items()}

        variant = {
            'genes': [gene_ids_by_name[line['Gene']]],
            'chrom': line['Chr'],
            'pos': int(line['Start'].replace(',', '')),
            'rsid': '.',
            'ref': line['Ref'],
            'alt': line['Var'],
            'eff_annotations' : [{
                'Gene_Name' : line['Gene'],
                'Amino_Acid_Change' : '/'.join(
                    line[nom] for nom in ('Protein', 'cDNA') if line[nom]),
                'Effect' : types.get(line['Type']) or line['Type'],
                # 'Transcript_BioType' : 'protein_coding',
                # 'Gene_Coding' : 'CODING',
                # 'Transcript_ID' : 'ENST00000345108',
                # 'Amino_Acid_length' : '321',
                # 'Effect_Impact' : 'MODIFIER',
            }],
            'allele_count': int(line['All carrier count']),
            'allele_num': ALLELE_NUM,
            'hom_count': int(line['Homoz']),
            'pop_acs': {
                'African American': int(line['AA']),
                'European American': int(line['EA']),
            },
            'pop_ans': POP_NUM,
            'quality_metrics': {'MQ': float(line['Max_MQ'])},
            'site_quality': float(line['Max_qual']),
            'transcripts': [line['Isoform']],
            'filter': 'PASS',
        }

        if line['Gene'] == 'RECQL':
            variant.update({
                'allele_num': ALLELE_NUM_RECQL,
                'pop_ans': POP_NUM_RECQL,
            })


        xpos = get_xpos(variant['chrom'], variant['pos'])
        variant.update({
            'variant_id': '-'.join((variant['chrom'], str(variant['pos']),
                                    variant['ref'], variant['alt'])),
            'allele_freq': float(variant['allele_count']) / variant['allele_num'],
            'xpos': xpos,
            'xstart': xpos,
            'xstop': xpos + len(variant['alt']) - len(variant['ref']),
        })

        yield variant


def get_mnp_data(mnp_file):
    header = mnp_file.readline().strip().split('\t')
    for line in mnp_file:
        data = dict(zip(header, line.split('\t')))
        if any(map(lambda x: x == 'True', data['QUESTIONABLE_PHASING'])): continue
        chroms = data['CHROM'].split(',')
        chrom = chroms[0]
        sites = data['SITES'].split(',')
        refs = data['REF'].split(',')
        alts = data['ALT'].split(',')
        for i, site in enumerate(sites):
            all_sites = zip(chroms, sites, refs, alts)
            all_sites.remove(all_sites[i])
            mnp = {}
            mnp['xpos'] = get_xpos(chrom, site)
            mnp['ref'] = refs[i]
            mnp['alt'] = alts[i]
            mnp['site2'] = '-'.join(all_sites[0])
            if len(all_sites) > 1:
                mnp['site3'] = all_sites[1]
            mnp['combined_codon_change'] = data['COMBINED_CODON_CHANGE']
            mnp['category'] = data['CATEGORY']
            mnp['number_samples'] = data['NSAMPS']
            yield mnp


def get_constraint_information(constraint_file):
    _, _, _, header = constraint_file.readline().strip().split(None, 3)
    header = header.split()
    for line in constraint_file:
        transcript, gene, chrom, info = line.strip().split(None, 3)
        transcript_info = dict(zip(header, map(float, info.split())))
        transcript_info['transcript'] = transcript.split('.')[0]
        yield transcript_info


def get_canonical_transcripts(canonical_transcript_file):
    for line in canonical_transcript_file:
        gene, transcript = line.strip().split()
        yield gene, transcript


def get_omim_associations(omim_file):
    for line in omim_file:
        fields = line.strip().split('\t')
        if len(fields) == 4:
            yield fields
        else:
            yield None


def get_genes_from_gencode_gtf(gtf_file):
    """
    Parse gencode GTF file;
    Returns iter of gene dicts
    """
    for line in gtf_file:
        if line.startswith('#'):
            continue
        fields = line.strip('\n').split('\t')

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
            'gene_name_upper': info['gene_name'].upper(),
            'chrom': chrom,
            'start': start,
            'stop': stop,
            'strand': fields[6],
            'xstart': get_xpos(chrom, start),
            'xstop': get_xpos(chrom, stop),
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
            'strand': fields[6],
            'xstart': get_xpos(chrom, start),
            'xstop': get_xpos(chrom, stop),
        }
        yield gene


def get_exons_from_gencode_gtf(gtf_file):
    """
    Parse gencode GTF file;
    Returns iter of transcript dicts
    """
    for line in gtf_file:
        if line.startswith('#'):
            continue
        fields = line.strip('\n').split('\t')

        if fields[2] not in ['exon', 'CDS', 'UTR']:
            continue

        chrom = fields[0][3:]
        feature_type = fields[2]
        start = int(fields[3]) + 1  # bed files are 0-indexed
        stop = int(fields[4]) + 1
        info = dict(x.strip().split() for x in fields[8].split(';') if x != '')
        info = {k: v.strip('"') for k, v in info.items()}
        transcript_id = info['transcript_id'].split('.')[0]
        gene_id = info['gene_id'].split('.')[0]

        exon = {
            'feature_type': feature_type,
            'transcript_id': transcript_id,
            'gene_id': gene_id,
            'chrom': chrom,
            'start': start,
            'stop': stop,
            'strand': fields[6],
            'xstart': get_xpos(chrom, start),
            'xstop': get_xpos(chrom, stop),
        }
        yield exon


def get_cnvs_from_txt(cnv_txt_file):
    """                                                                                                                                                                                                                                  
    Parse gencode txt file;                                                                                                                                                                                                              
    Returns iter of gene dicts                                                                                                                                                                                                           
    """
    header = cnv_txt_file.next() # gets rid of the header                                                                                                                                                                                
    #print header                                                                                                                                                                                                                        
    for line in cnv_txt_file:

        fields = line.rsplit()
        transcript = fields[0]
        gene = fields[1]
        chrom = fields[2]
        start = int(fields[3])
        stop = int(fields[4])
        del0 = int(fields[5])
        del60 = int(fields[6])
        dup0 = int(fields[7])
        dup60 = int(fields[8])
        delpop0 = fields[9]
        delpop60 = fields[10]
        duppop0 = fields[11]
        duppop60 = fields[12]
        

        #find gene from DB.genes, get ID                                                                                                                                                                                                 
        #find exon of that gene that this CNV referes to from db.exons, get ID                                                                                                                                                           
        #add the object reference to the cnv dict.                                                                                                                                                                                       
        cnv = {
            'transcript': transcript,
            'gene': gene,
            'chrom': chrom,
            'start': start,
            'stop': stop,
            'del0': del0,
            'dup0': dup0,
            'dup60': dup60,
            'del60' : del60,
            'delpop0' : delpop0,
            'delpop60' : delpop60,
            'duppop0' : duppop0,
            'duppop60' : duppop60,
            'xstart': get_xpos(chrom, start),
            'xstop': get_xpos(chrom, stop),
        }
        yield cnv


def get_cnvs_per_gene(cnv_gene_file):
    header = cnv_gene_file.next() # gets rid of the header                                                                                                                                                                               
    for line in cnv_gene_file:

        fields = line.rsplit()
        gene = fields[0]
        symbol = fields[1]
        del0 = int(fields[2])
        dup0 = int(fields[3])
        cnv0 = int(fields[4])
        del60 = int(fields[5])
        dup60 = int(fields[6])
        cnv60 = int(fields[7])
        del_score = float(fields[8])
        dup_score = float(fields[9])
        cnv_score = float(fields[10])
        rank = int(fields[11])

        cnv_gene = {
            'gene': gene,
            'symbol': symbol,
            'del0': del0,
            'dup0': dup0,
            'cnv0': cnv0,
            'del60': del60,
            'dup60': dup60,
            'cnv60' : cnv60,
            'del_score': del_score,
            'dup_score': dup_score,
            'cnv_score': cnv_score,
            'rank': rank,
            }
        yield cnv_gene




def get_dbnsfp_info(dbnsfp_file):
    """
    Parse dbNSFP_gene file;
    Returns iter of transcript dicts
    """
    header = dbnsfp_file.next().split('\t')
    fields = dict(zip(header, range(len(header))))
    for line in dbnsfp_file:
        line = line.split('\t')
        other_names = line[fields["Gene_old_names"]].split(';') if line[fields["Gene_old_names"]] != '.' else []
        if line[fields["Gene_other_names"]] != '.':
            other_names.extend(line[fields["Gene_other_names"]].split(';'))
        gene_info = {
            'gene_name': line[fields["Gene_name"]],
            'ensembl_gene': line[fields["Ensembl_gene"]],
            'gene_full_name': line[fields["Gene_full_name"]],
            'gene_other_names': other_names
        }
        yield gene_info


def get_snp_from_dbsnp_file(dbsnp_file):
    for line in dbsnp_file:
        fields = line.split('\t')
        if len(fields) < 3: continue
        rsid = int(fields[0])
        chrom = fields[1].rstrip('T')
        if chrom == 'PAR': continue
        start_str = fields[2]
        # skip initial rows like:
        # 869027110	1		0			
        if not start_str: continue
        start = int(start_str) + 1
        snp = {
            'xpos': get_xpos(chrom, start),
            'rsid': rsid
        }
        yield snp
