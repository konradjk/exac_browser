# Adapted from Daniel Birnbaum's histogram script

import argparse
import gzip
import pipes
import sys
from collections import Counter
import numpy

metrics = ['DP', 'GQ']

def main(args):
    f = gzip.open(args.vcf) if args.vcf.endswith('.gz') else open(args.vcf)

    if args.output is None: args.output = args.vcf.replace('.vcf', '.hist.vcf')
    if not args.output.endswith('.gz'): args.output += '.gz'

    pipe = pipes.Template()
    pipe.append('bgzip -c /dev/stdin', '--')
    g = pipe.open(args.output, 'w')

    header = None
    for line in f:
        line = line.strip()

        # Reading header lines to get VEP and individual arrays
        if line.startswith('#'):
            line = line.lstrip('#')
            if line.startswith('CHROM'):
                header = line.split()
                header = dict(zip(header, range(len(header))))
            continue

        if header is None:
            print >> sys.stderr, "VCF file does not have a header line (CHROM POS etc.). Exiting."
            sys.exit(1)

        fields = line.split('\t')
        # Pull out annotation info from INFO and ALT fields
        new_info = fields[header['INFO']].rstrip(';')

        for metric in metrics:
            data = get_histogram_for_variant(line, metric)
            midpoints, hist = data
            new_info += ';%s_MID=' % (metric) + '|'.join(map(str, midpoints))
            new_info += ';%s_HIST=' % (metric) + '|'.join(map(str, hist))

        fields[header['INFO']] = new_info
        g.write('\t'.join(fields) + '\n')

    f.close()
    g.close()


def convert_to_int(val):
    """
    Converts string to int if possible, otherwise returns initial string
    """
    try:
        return int(val)
    except ValueError:
        pass
    try:
        return float(val)
    except ValueError:
        return val


def get_histogram_for_variant(vcf_line, metric="DP", num_bins=40, midpoints=True, variants_only=False):
    vcf_line = vcf_line.strip('\n')
    if vcf_line.startswith('#'):
        return None
    else:
        fields = vcf_line.split('\t')
        # alts = fields[4].split(',')
        try:
            idx = fields[8].split(':').index(metric)
        except Exception, e:
            return None

        distr = []
        # get distribution for metric
        for sample in fields[9:]:
            # This is only DP/GQ for now
            sample_info = sample.split(':')
            if idx < len(sample_info) and sample_info[idx] != '.':
                distr.append(sample_info[idx])

        mids, hist = get_hist_from_distribution(distr, midpoints, num_bins)
        return map(str, mids), map(str, hist)


def get_hist_from_distribution(distr, midpoints, num_bins):
    distr = [convert_to_int(x) for x in distr]
    if any([type(x) == str for x in distr]):
        c = Counter(distr)
        counts = zip(*c.items())
        return counts
    else:
        hist = numpy.histogram(distr, bins=num_bins)
        if midpoints:
            edges = hist[1]
            mids = [(edges[i]+edges[i+1])/2 for i in range(len(edges)-1)]
            return mids, hist[0]
        else:
            return hist[1], hist[0]

if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('--vcf', '--input', '-i', help='Input VCF file; may be gzipped', required=True)
    parser.add_argument('--output', '-o', help='Output VCF file; may be gzipped')
    args = parser.parse_args()
    main(args)