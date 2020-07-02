"""Contains data IO and processing functionalities and kernel functions specific to Master thesis of Pia Rautenstrauch.
"""

# imports
import logging

import numpy as np
import pandas as pd
import pyreadr

from seak.data_loaders import CovariatesLoader

# set logging configs
logging.basicConfig(format='%(asctime)s - %(lineno)d - %(message)s')


class StefanReplicateCovariatesLoader(CovariatesLoader):
    __slots__ = ['phenotype_of_interest']

    def __init__(self, phenotype_of_interest):
        """Needs to be built once per set and phenotype to test."""
        path_to_covariates = "/mnt/30T/ukbiobank/derived/projects/kernels_VEP/ukb27451_obesity_diabetes_final.csv"
        self.phenotype_of_interest = phenotype_of_interest
        self.cov = pd.read_csv(path_to_covariates, index_col=0)
        self.cov['cov_nindex'] = self.cov.index
        self.cov = self.cov.rename(columns={'eid': 'iid'})

        # Index of related individuals provided by UK Biobank from Stefan --> we want to exclude those
        # Load RData file
        eid_unrel = pyreadr.read_r("/mnt/30T/ukbiobank/derived/projects/kernels_VEP/idx_unrel_exome_seq.RData")
        eid_unrel = eid_unrel['idx_unrel']
        eid_unrel = eid_unrel['idx_unrel'].tolist()

        # Only consider individuals for which information on all covariates and phenotype is available!
        # As Stefan did in the script: UKBiobank_2_SKAT_SKATO
        # The sample index is called eid instead of iid here!

        # Drop all individuals with incomplete covariate or phenotype information
        self.cov = self.cov[['iid', 'cov_nindex', self.phenotype_of_interest, 'age', 'sex', 'smoke', 'income',
                             'Caucasian', 'genet_PC_1', 'genet_PC_2', 'genet_PC_3', 'genet_PC_4', 'genet_PC_5',
                             'genet_PC_6', 'genet_PC_7', 'genet_PC_8', 'genet_PC_9', 'genet_PC_10']].dropna()

        # Remove strongly related individuals (defined by Stefans R script, I only load the index)
        self.cov = self.cov[~(self.cov['iid'].isin(eid_unrel))]
        print('# of samples with complete genotype and covariate information: ', self.cov.shape[0])
        self.cov['iid'] = self.cov['iid'].astype(str)
        self.cov.set_index(keys='iid', inplace=True, drop=False)

    def update_individuals(self, iids):
        """Set individuals to include.

        This also changes the order to match the order in :attr:`iids`.

        :param pandas.Index iids: Individual ids to include
        """
        self.cov = self.cov.loc[iids]

    def get_one_hot_covariates_and_phenotype(self, test_type):
        """Returns numpy.ndarray of the phenotype and one hot encoded covariates with invariant covariates
        removed and a bias column.

        Make sure :func:`update_cov` was called beforehand, such that :attr:`self.cov` has the same order as the genotypes.
        Internally calls function :func:`pandas.get_dummies` which converts categorical variable into indicator variables.
        Should be repeated per strand.

        :return: data for phenotype and one hot encoded covariates with invariant covariates removed and bias column
        :rtype: numpy.ndarray
        """
        print('Get one hot covariates and phenotype')
        one_hot_covariates = pd.get_dummies(self.cov[['age', 'sex', 'smoke', 'income', 'genet_PC_1', 'genet_PC_2',
                                                      'genet_PC_3', 'genet_PC_4', 'genet_PC_5', 'genet_PC_6',
                                                      'genet_PC_7', 'genet_PC_8', 'genet_PC_9', 'genet_PC_10']],
                                            prefix_sep='_', drop_first=True)

        # Drop invariant covariates
        one_hot_covariates = one_hot_covariates.loc[:, one_hot_covariates.apply(pd.Series.nunique) != 1]
        one_hot_covariates = np.asarray(one_hot_covariates)
        # Add offset/bias column
        X = np.hstack((np.ones((one_hot_covariates.shape[0], 1)), one_hot_covariates))
        if test_type == 'logit':
            phenotype = np.asarray(pd.get_dummies(self.cov[self.phenotype_of_interest], dtype=float, drop_first=True))
            if phenotype.shape[1] != 1:
                logging.ERROR('It seems like your phenotype was not encoded binary.')
            phenotype = phenotype.reshape(len(phenotype), 1)
        else:
            phenotype = np.asarray(self.cov[self.phenotype_of_interest])
            phenotype = phenotype.reshape(len(phenotype), 1)
        return phenotype, X

    def get_iids(self):
        """Returns all individual ids.

        :return:
        :rtype: pandas.Index
        """
        return self.cov.index


class StefanReplicateCovariatesLoader_with_stratification(StefanReplicateCovariatesLoader):
    """Loads all UK biobank exome samples except related individuals."""
    def __init__(self, phenotype_of_interest):
        """Needs to be built once per set and phenotype to test."""
        path_to_covariates = "/mnt/30T/ukbiobank/derived/projects/kernels_VEP/ukb27451_obesity_diabetes_final.csv"
        self.phenotype_of_interest = phenotype_of_interest
        self.cov = pd.read_csv(path_to_covariates, index_col=0)
        self.cov['cov_nindex'] = self.cov.index
        self.cov = self.cov.rename(columns={'eid': 'iid'})

        # Index of related individuals provided by UK Biobank from Stefan --> we want to exclude those
        # Load RData file
        eid_unrel = pyreadr.read_r("/mnt/30T/ukbiobank/derived/projects/kernels_VEP/idx_unrel_exome_seq.RData")
        eid_unrel = eid_unrel['idx_unrel']
        eid_unrel = eid_unrel['idx_unrel'].tolist()

        # Only consider individuals for which information on all covariates and phenotype is available!
        # As Stefan did in the script: UKBiobank_2_SKAT_SKATO

        # Drop all individuals with incomplete covariate or phenotype information
        self.cov = self.cov[['iid', 'cov_nindex', self.phenotype_of_interest, 'age', 'sex', 'smoke', 'income',
                             'genet_PC_1', 'genet_PC_2', 'genet_PC_3', 'genet_PC_4', 'genet_PC_5',
                             'genet_PC_6', 'genet_PC_7', 'genet_PC_8', 'genet_PC_9', 'genet_PC_10']].dropna()

        self.cov = self.cov[~(self.cov['iid'].isin(eid_unrel))]
        print('# of samples with complete genotype and covariate information: ', self.cov.shape[0])
        self.cov['iid'] = self.cov['iid'].astype(str)
        self.cov.set_index(keys='iid', inplace=True, drop=False)


class StefanReplicateCovariatesLoader_with_stratification_and_relatedness(StefanReplicateCovariatesLoader):
    """Loads all UK biobank exome samples with complete information."""
    def __init__(self, phenotype_of_interest):
        """Needs to be built once per set and phenotype to test."""
        path_to_covariates = "/mnt/30T/ukbiobank/derived/projects/kernels_VEP/ukb27451_obesity_diabetes_final.csv"
        self.phenotype_of_interest = phenotype_of_interest
        self.cov = pd.read_csv(path_to_covariates, index_col=0)
        self.cov['cov_nindex'] = self.cov.index
        self.cov = self.cov.rename(columns={'eid': 'iid'})

        # Only consider individuals for which information on all covariates and phenotype is available!
        # As Stefan did in the script: UKBiobank_2_SKAT_SKATO

        # Drop all individuals with incomplete covariate or phenotype information
        self.cov = self.cov[['iid', 'cov_nindex', self.phenotype_of_interest, 'age', 'sex', 'smoke', 'income',
                             'genet_PC_1', 'genet_PC_2', 'genet_PC_3', 'genet_PC_4', 'genet_PC_5',
                             'genet_PC_6', 'genet_PC_7', 'genet_PC_8', 'genet_PC_9', 'genet_PC_10']].dropna()

        print('# of samples with complete genotype and covariate information: ', self.cov.shape[0])
        self.cov['iid'] = self.cov['iid'].astype(str)
        self.cov.set_index(keys='iid', inplace=True, drop=False)


# Kernel functions
def diffscore_max_unscaled(G, V):
    """Uses the largest absolute variant effect predictions (veps) per SNV as linear weights.

    Linear weighted kernel

    :param numpy.ndarray G: SNVs to be tested (genotypes), :math:`nxm` (:math:`n:=` number of individuals, :math:`m:=` number of SNVs)
    :param numpy.ndarray V: veps :math:`mxk` (:math:`m:=` number of SNVs, :math:`k:=` number of veps); dimension :math:`m` needs to be equal in G and V
    :return: linear weighted genotype matrix (:math:`GxV`)
    :rtype: numpy.ndarray
    """
    V = np.abs(V)
    # Maximum out of any diffscore prediction
    V = np.amax(V, axis=1)
    V = np.diag(V, k=0)
    return np.matmul(G, V)


def diffscore_max_collapsed(G, V):
    """Uses the largest absolute variant effect predictions (veps) per SNV as linear weights.

    Linear weighted kernel

    :param numpy.ndarray G: SNVs to be tested (genotypes), :math:`nxm` (:math:`n:=` number of individuals, :math:`m:=` number of SNVs)
    :param numpy.ndarray V: veps :math:`mxk` (:math:`m:=` number of SNVs, :math:`k:=` number of veps); dimension :math:`m` needs to be equal in G and V
    :return: linear weighted genotype matrix (:math:`GxV`)
    :rtype: numpy.ndarray
    """
    V = np.abs(V)
    # Maximum out of any diffscore prediction
    V = np.amax(V, axis=1)
    return np.matmul(G, V).reshape(G.shape[0], 1)


def diffscore_max_unscaled_sqrt(G, V):
    """Uses the largest absolute variant effect predictions (veps) per SNV as linear weights.

    Linear weighted kernel

    :param numpy.ndarray G: SNVs to be tested (genotypes), :math:`nxm` (:math:`n:=` number of individuals, :math:`m:=` number of SNVs)
    :param numpy.ndarray V: veps :math:`mxk` (:math:`m:=` number of SNVs, :math:`k:=` number of veps); dimension :math:`m` needs to be equal in G and V
    Dimension m needs to be equal in G and V
    :return: linear weighted genotype matrix (:math:`GxV`)
    :rtype: numpy.ndarray
    """
    V = np.sqrt(np.abs(V))
    # Maximum out of any diffscore prediction
    V = np.amax(V, axis=1)
    V = np.diag(V, k=0)
    return np.matmul(G, V)


def diffscore_max_collapsed_sqrt(G, V):
    """Uses the largest absolute variant effect predictions (veps) per SNV as linear weights.

    Linear weighted kernel

    :param numpy.ndarray G: SNVs to be tested (genotypes), :math:`nxm` (:math:`n:=` number of individuals, :math:`m:=` number of SNVs)
    :param numpy.ndarray V: veps :math:`mxk` (:math:`m:=` number of SNVs, :math:`k:=` number of veps); dimension :math:`m` needs to be equal in G and V
    :return: linear weighted genotype matrix (:math:`GxV`)
    :rtype: numpy.ndarray
    """
    V = np.sqrt(np.abs(V))
    # Maximum out of any diffscore prediction
    V = np.amax(V, axis=1)
    return np.matmul(G, V).reshape(G.shape[0], 1)




