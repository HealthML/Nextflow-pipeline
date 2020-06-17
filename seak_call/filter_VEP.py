import numpy as np
import pandas as pd
import os
import argparse

    # awk '$NF ~ /IMPACT=HIGH/' /home/Aliki.Zavaropoulou/UKbiobank/derived/projects/kernels_VEP/vep_SPB_out/vep_ensembl/v2/output/vep/ukb_SPB.vep.vcf > /home/Aliki.Zavaropoulou/UKbiobank/derived/projects/LMM-Lasso/deepWAS/ukb_SPB_filteredvariants.vep.vcf
def main(input_file, output_file):
    # positions in vcf files are 1-based in "UCSC BED like format" they are 0-based and pre-filtered variants
    # in_path = '/home/Aliki.Zavaropoulou/UKbiobank/derived/projects/kernels_VEP/vep_SPB_out/vep_ensembl/v1/output/my_output.vcf' # note: it's not actually a VCF file
    # out_path = './ensembl_vep_necessary'
    # os.makedirs(out_path, exist_ok=True)
    # read the VEP output
    #colnames = ['Uploaded_variation', 'Location', 'Allele', 'Gene', 'Feature', 'Feature_type', 'Consequence', 'cDNA_position', 'CDS_position', 'Protein_position',  'Amino_acids',  'Codons',  'Existing_variation', 'Extra' ]
    colnames = ['Uploaded_variation', 'Location', 'Allele', 'Gene', 'Consequence']
    cols = [0, 1, 2, 3, 6]
    df = pd.read_csv(input_file, sep='\t', comment='#', header=None, usecols=cols, names=colnames)

    # LOF categories; We use the same as in Cirulli et. al, however we don't have any frameshift variants because in v1 we filtered out indels
    lof_anno = pd.Series(['stop_lost','start_lost','splice_donor_variant','frameshift_variant','splice_acceptor_variant','stop_gained'])

    is_lof = df['Consequence'].str.split(',', expand=True)
    is_lof = np.any(is_lof.apply(lambda x: x.isin(lof_anno), axis=0).values, axis=1)
    df = df[is_lof]

    df_summary = df.groupby(['Uploaded_variation','Location','Allele','Gene']).size()
    df_summary = df_summary.reset_index()
    df_summary.rename(columns={0:'n_affected'}, inplace=True)
    # print(df_summary.shape)   ->  (264528, 5)
    # print(df_summary.head())
    df_summary.to_csv(output_file, sep='\t', index=False)


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("-i", help="input file name", default = '/home/Aliki.Zavaropoulou/UKbiobank/derived/projects/kernels_VEP/test_pipeline/ukb_SPB_filteredvariants.vep.vcf', type = str)
    parser.add_argument("-o", help="output file name", default = '/home/Aliki.Zavaropoulou/pipeline/Nextflow-pipeline/seak_call/ensembl_vep_necessary/LOF_filtered.tsv', type = str)
    args = parser.parse_args()

    # Call main routine
    main(args.i, args.o)
