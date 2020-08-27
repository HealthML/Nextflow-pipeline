# comments and maybe modifications by Aliki Zavaropoulou, June 09
# ----------------------------------------------------------------------------------------------------------------------
# May 18, 2020
# Command line interface to loop over chromosomes
# Include DeepRiPE predictions for specific RBPs, restrict analysis to the corresponding variants.
# Condition on protein LOF mutations
# use these same variants with without weights (linear kernel)
# use local collapsing kernel
# Use of exome data set excluding related individuals
# Test: noK
# ----------------------------------------------------------------------------------------------------------------------

import time
import pandas as pd
import logging
import numpy as np
import argparse
import os
import sys
#import src.seak.data_loaders
#temporary way to load the seak library until it is available in pip
# sys.path.insert(0, '/home/Aliki.Zavaropoulou/pipeline/Nextflow-pipeline/seak/src/seak')
# import data_loaders
# import kernels
# import scoretest

from seak import data_loaders
from seak import kernels
from seak import scoretest

from argparse import ArgumentParser
from numpy.linalg import LinAlgError

# set logging configs
logging.basicConfig(format='%(asctime)s - %(lineno)d - %(message)s')


def get_args():

    parser = ArgumentParser()

    parser.add_argument('-pheno', '--phenotype', type=str, required=True,
                       choices={'ApoA', 'ApoB', 'IGF1', 'CRP'})
    parser.add_argument('-i', '--input', default = 'LOF_filtered.tsv', type = str)
    parser.add_argument('-ref', '--referencegenome', default = '/mnt/dsets/reference_genomes/ensembl/Homo_sapiens.GRCh38.genes.bed', type = str)
    parser.add_argument('-ukbdir', '--ukbiobankdirectory', default = '/home/Aliki.Zavaropoulou/UKbiobank', type = str)
    parser.add_argument('-g', '--geno', default = '${/home/Aliki.Zavaropoulou/UKbiobank/derived/projects/kernels_VEP/test_pipeline}/ukb_SPB_50k_exome_seq_filtered', type = str)
    parser.add_argument('-l', '--lofpath', default = '/home/Aliki.Zavaropoulou/UKbiobank/derived/projects/kernels_VEP/test_pipeline/seaktsv/')
    parser.add_argument('-o', '--outputpath', default = '/home/Aliki.Zavaropoulou/UKbiobank/derived/projects/kernels_VEP/test_pipeline/seaktsv/')
    parser.add_argument('-cov', '--covariatespath', default = '/home/Aliki.Zavaropoulou/UKbiobank/derived/projects/kernels_VEP/master_thesis_pia/blood_biochemistry_phenotypes_and_covariates/blood_biochemistry_pheno_and_cov_unrelated_gaussian_quantile_transformed.csv' )
    # parser.add_argument('-chr', '--chromosome', type=str, required=True)
    # parser.add_argument('-minvep', default=0.1, type=float)
    # parser.add_argument('-maxdist', default=150, type=int)

    return parser.parse_args()

def weights_sign_from_V(V):
    if V.shape[1] == 1:
        return np.sqrt(np.abs(V[:,0])), np.sign(V[:,0])
    else:
        V = V[np.arange(V.shape[0]), np.argmax(np.abs(V), axis=1)]
        return np.sqrt(np.abs(V)), np.sign(V)

def main():

    args = get_args()

    # parameters
    ukb_dir = args.ukbiobankdirectory
    test_type = 'noK'
    drop_non_numeric_chromosomes = True
    max_maf = 0.001
    path_to_regions_UCSC_BED = args.referencegenome
    pheno_short = args.phenotype

    short_to_pheno = {'ApoA':'Gaussian_quantile_transform(Apolipoprotein A)', 'ApoB':'Gaussian_quantile_transform(Apolipoprotein B)', 'IGF1': 'Gaussian_quantile_transform(IGF-1)', 'CRP':'Gaussian_quantile_transform(C-reactive protein)'}
    short_to_dir = {'ApoA': 'ApoA/', 'ApoB': 'ApoB/', 'IGF1': 'IGF1/', 'CRP': 'CRP/'}

    phenotype_of_interest = short_to_pheno[pheno_short]
    experiment = 'noK_unrelated_filtered_protLOF_maf0.001_burden'

    # protein LOF variants - here is probably the connection to the VEP procedure (initially: /home/Aliki.Zavaropoulou/UKbiobank/derived/projects/kernels_VEP/vep_SPB_out/v2/ensembl_vep/lof_filtered.tsv)
    ensemblvep_df = pd.read_csv(args.lofpath + args.input , sep='\t')
    #ensemblvep_df = pd.read_csv(ukb_dir + '/derived/projects/kernels_VEP/test_pipeline/seaktsv/' + args.input , sep='\t')
    # ensemblvep_df = pd.read_csv(ukb_dir + '/derived/projects/kernels_VEP/vep_SPB_out/v2/ensembl_vep/' + args.input , sep='\t')

    # SNPs we want to test using the whole genome .bed file generated from process pling_1
    path_to_plink_geno_from_pipeline = args.geno    #not used
    # this variable is defined twice and here it is probably not used
    path_to_plink_geno_files_with_prefix = ukb_dir+ '/derived/projects/kernels_VEP/test_pipeline/ukb_FE_50k_exome_seq_filtered/ukb_FE_50k_exome_seq_filtered'
    plinkloader_G1 = data_loaders.VariantLoaderSnpReader(path_or_bed=path_to_plink_geno_files_with_prefix+'.bed')

    for chromosome in range(1, 23):

        print('Experiment, chromosome:', experiment, str(chromosome))

        chrom_to_load = str(chromosome)

        # TODO: change the output directory ; changed from './results/experiments/blood_biochemistry/' -> './results/' -> '/home/Aliki.Zavaropoulou/UKbiobank/derived/projects/kernels_VEP/test_pipeline/results/'
        output_prefix = args.outputpath + experiment + '/' + short_to_dir[pheno_short] + '/' + chrom_to_load
        os.makedirs(os.path.dirname(output_prefix), exist_ok=True)

        # plink file prefix - this variable is defined twice
        path_to_plink_geno_files_with_prefix = ukb_dir+ '/derived/projects/kernels_VEP/ukb_SPB_50k_exome_seq_filtered_chr' + chrom_to_load

        # protein LOF variants
        ensemblvep_loader = data_loaders.EnsemblVEPLoader(ensemblvep_df['Uploaded_variation'], ensemblvep_df['Location'], ensemblvep_df['Gene'], data=None)

        path_to_covariates = args.covariatespath
        covariate_column_names = ['Age at recruitment', 'Sex', 'BMI', 'Smoking status', 'PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6', 'PC7', 'PC8', 'PC9', 'PC10']

        # SNPs we want to test - TODO: path_to_plink_ ... should be an input in the procedure and in arg.parser
        plinkloader_G1 = data_loaders.VariantLoaderSnpReader(path_or_bed=path_to_plink_geno_files_with_prefix+'.bed')

        # covariates
        covariate_loader_csv = data_loaders.CovariatesLoaderCSV(phenotype_of_interest=phenotype_of_interest, path_to_covariates=path_to_covariates, covariate_column_names=covariate_column_names)

        data_loaders.intersect_and_update_datasets(test=test_type, variantloader=plinkloader_G1, covariateloader=covariate_loader_csv, annotationloader=ensemblvep_loader)

        print('{} protein LOF variants left after filtering'.format(len(ensemblvep_loader.get_vids())))

        Y, X = covariate_loader_csv.get_one_hot_covariates_and_phenotype(test_type=test_type)

        # initialize null model
        null_model = scoretest.ScoretestNoK(Y, X)

        # get regions, set up iterator to loop over
        #path_to_regions_UCSC_BED = '/mnt/dsets/reference_genomes/ensembl/Homo_sapiens.GRCh38.genes.bed'
        ucsc_region_loader = data_loaders.BEDRegionLoader(path_to_regions_UCSC_BED=path_to_regions_UCSC_BED, chrom_to_load=chrom_to_load, drop_non_numeric_chromosomes=drop_non_numeric_chromosomes)


        def load_and_proc_geno(interval):
            # load temp_G1
            temp_genotypes, temp_vids, temp_pos = plinkloader_G1.genotypes_by_region(interval, return_pos=True)

            try:
                V1 = ensemblvep_loader.anno_by_interval(interval, gene=interval['name'].split('_')[0])
            except KeyError:
                return (None, None)

            if V1.index.empty:
                return (None, None)

            vids = V1.index.get_level_values('vid')
            temp_genotypes, temp_vids = plinkloader_G1.genotypes_by_id(vids, return_pos=False)

            G1, vids = plinkloader_G1.preprocess_genotypes(temp_genotypes, temp_vids, recode_maf=False, invert_encoding=True, impute_mean=False, max_maf = max_maf)

            if G1 is None:
                return (None, None)

            G1_burden = (np.nansum(G1, axis=1, keepdims=True) > 0.).astype(float) # indicator: 1 if at least one variant, 0 otherwise

            return (G1_burden, vids)


        def testGV(GV):
            p_uncond = null_model.pv_alt_model(GV)
            return p_uncond


        def test(interval):

            (G1, vids) = load_and_proc_geno(interval)

            pval_dict = {}

            if G1 is None:
                return None, None

            # linear kernel (no weights)
            pval_dict['linb'] = testGV(G1)

            pval_dict['nSNP'] = len(vids)

            return pval_dict, vids

        results_out = output_prefix + '_results.tsv'
        ids_path = os.path.dirname(output_prefix)+'ids/'
        os.makedirs(ids_path, exist_ok=True)

        with open(results_out, 'w') as out:
            out.write('gene\tlinb\tnSNP\n')

        for region in ucsc_region_loader:

            pv, vids = test(region)

            if pv is None:
                continue

            pvdf = pd.DataFrame(pv, index=[region['name']])
            pvdf.to_csv(sep='\t',index_label='gene', path_or_buf=results_out, header=False, mode='a')

            if pvdf['linb'].values <= 0.005:
                # export variant IDs for most significant hits
                with open(ids_path + region['name']+'_ids.txt', 'w') as out:
                    out.write(','.join(vids)+'\n')


if __name__ == '__main__':
    main()
