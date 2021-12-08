import dms_variants.codonvarianttable
import dms_variants.illuminabarcodeparser
import argparse
import sys, os
import re
import pandas as pd
import Bio.SeqIO

illu_params = {
    'bclen': 26, 
    'upstream': '', 
    'downstream': '', 
    'minq': 20, 
    'upstream_mismatch': 0,
    'downstream_mismatch': 0}

parser = argparse.ArgumentParser(description="Step 1: count barcodes for fastq file.")

parser.add_argument('--input_fastq', '-i', type=str)
parser.add_argument('--output_dir', '-o', type=str)
# parser.add_argument('--library', '-l', type=str)
parser.add_argument('--table', '-t', type=str)
parser.add_argument('--wildtype', '-w', type=str)

def build_variant_table(table_path, wt_seq_path):
    wt_seqrecord = Bio.SeqIO.read(wt_seq_path, 'fasta')
    geneseq = str(wt_seqrecord.seq)
    primary_target = wt_seqrecord.name
    variants = dms_variants.codonvarianttable.CodonVariantTable(
               geneseq = geneseq,
               barcode_variant_file = table_path,
               substitutions_are_codon = True,
               substitutions_col = 'codon_substitutions',
               primary_target = primary_target)
    return variants, primary_target

def count_variants(variants, lib, sample, fastq_data, illumina_barcode_parser_params, output_dir):
    # variants: returned object of dms_variants.codonvarianttable
    parser = dms_variants.illuminabarcodeparser.IlluminaBarcodeParser(
             valid_barcodes = variants.valid_barcodes(lib),
             **illumina_barcode_parser_params)
    counts, fates = parser.parse(fastq_data, add_cols = {'library':lib, 'sample': sample})

    counts_file = os.path.join(output_dir, '_'.join([sample,'counts'])+'.csv')
    fates_file = os.path.join(output_dir, '_'.join([sample,'fates'])+'.csv')

    counts.to_csv(counts_file)
    fates.to_csv(fates_file)

    return counts, fates

def main():
    args = parser.parse_args()
    sample_name = args.input_fastq.split('/')[-1].split('.')[0]
    library = sample_name.split('_')[-1]
    #library = re.search(r"lib[1-2]", sample_name).group()
    variants, primary_target = build_variant_table(args.table, args.wildtype)
    output_dir = os.path.join(args.output_dir, sample_name)
    os.makedirs(output_dir, exist_ok=True)

    sys.stdout.write('Start parsing barcodes...\n')
    if os.path.exists(os.path.join(output_dir, sample_name+'_counts.csv')) and os.path.exists(os.path.join(output_dir, sample_name+'_fates.csv')):
        counts = pd.read_csv(os.path.join(output_dir, sample_name+'_counts.csv'), index_col=0)
        sys.stdout.write('Skipped! '+sample_name+'\n')
    else:
        counts, fates = count_variants(variants, library, sample_name, fastq_data = args.input_fastq, illumina_barcode_parser_params=illu_params, output_dir=output_dir)
        sys.stdout.write('Done.\n')
    
    if not os.path.exists(os.path.join(output_dir, sample_name+'_variant_counts.csv')):
        variants.add_sample_counts_df(counts)
        variant_counts = variants.variant_count_df
        variant_counts.to_csv(os.path.join(output_dir, sample_name+'_variant_counts.csv'), index = False)

if __name__ == '__main__':
    main()
