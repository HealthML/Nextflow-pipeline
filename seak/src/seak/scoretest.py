"""Contains classes for score-based set-association tests.

A single kernel and two-kernel score-based test is available for a linear link function (continuous phenotypes) (:class:`ScoretestNoK` and :class:`Scoretest2K`).
With the second kernel :math:`K_0` correcting for population structure.
A single kernel score-based test is available for a logistic link function (binary phenotypes) (:class:`ScoretestLogit`)

Null model single kernel:

.. math:: y = {\\alpha} X + {\epsilon}

Alternative model single kernel:

.. math:: y = {\\alpha} X + {\\gamma} G_1 + {\\epsilon}

Null model two-kernel:

.. math:: y = {\\alpha} X + {\\beta} G_0 + {\\epsilon}

Alternative model two-kernel:

.. math:: y = {\\alpha} X + {\\beta} G_0 + {\\gamma} G_1 + {\\epsilon}

With: :math:`X`: covariates, dimension :math:`nxc` (:math:`n:=` number of individuals, :math:`c:=` number of covariates),
:math:`G_0`: SNVs to construct the background kernel :math:`K_0` from, correcting for population structure, dimensions :math:`nxm` (:math:`n:=` number of individuals, :math:`m:=` number of SNVs/variants),
:math:`G_1`: set of SNVs to test, dimensions :math:`nxk` (:math:`n:=` number of individuals, :math:`k:=` number of SNVs in set to test association for).

.. note:: For all classes of the module :mod:`scoretest` the class attributes are instance attributes! This is a bug in the automatic documentation.

.. note::
   The source code of the mathematical implementation is adopted from the `FaST-LMM <https://pypi.org/project/fastlmm/>`_ Python package  from Microsoft Corporation.
   Source code can be found in :mod:`fastlmm.association.score`.

"""

# imports
import logging
import time

import numpy as np
import scipy as sp
import scipy.linalg as LA
import statsmodels.api as sm

from seak.mingrid import minimize1D

# set logging configs
logging.basicConfig(format='%(asctime)s - %(lineno)d - %(message)s')


class Scoretest:
    """Superclass for all score-based set-association tests.

    Defines all common attributes and implements methods that all subclasses share, such as methods for p-value
    computation for alternative :func:`pv_alt_model` model.

    Interface for Scoretest subclasses.

    :ivar numpy.ndarray Y: :math:`nx1` matrix containing the phenotype (:math:`n:=` number of individuals)
    :ivar numpy.ndarray X: :math:`nxc` matrix containing the covariates (:math:`n:=` number of individuals, :math:`c:=` number of covariates)
    :ivar int N: number of individuals
    :ivar int P: number of phenotypes
    :ivar int D: number of covariates
    :ivar int Neff: effective sample size (number of individuals) as used in calculations of unbiased estimators
    """

    __slots__ = ["Y", "X", "N", "P", "D", "Neff"]

    def __init__(self, phenotypes, covariates):
        """Constructor.

        If no covariates are given or if the covariates do not contain a bias column, a bias column is added.
        Sets up all attributes common to all score-based set-association tests.

        :param numpy.ndarray phenotypes: :math:`nx1` matrix containing the phenotype (:math:`n:=` number of individuals)
        :param numpy.ndarray covariates: :math:`nxc` matrix containing the covariates (:math:`n:=` number of individuals, :math:`c:=` number of covariates)
        """
        self.Y = self._set_phenotypes(phenotypes)  # phenotypes
        self.N = self.Y.shape[0]  # number of individuals
        self.P = self.Y.shape[1]  # number of phenotypes
        self.X = self._set_covariates(covariates)  # covariates
        self.D = self.X.shape[1]  # num of covariates
        self.Neff = self.N - self.D  # unbiased estimator for variance

        if self.P != 1:
            logging.ERROR('More than one phenotype given.')

        if self.X.shape[0] != self.N:
            logging.ERROR('Number of individuals in phenotype and covariates does not match.')

    def _set_phenotypes(self, phenotypes):
        """Casts phenotypes to two dimensions."""
        if phenotypes.ndim == 1:
            phenotypes = phenotypes.reshape(len(phenotypes), 1)
        return phenotypes

    def _set_covariates(self, covariates):
        """Appends bias (offset) column to covariates if not present."""
        if covariates is None:
            X = np.ones((self.N, 1))
        elif not self.has_bias(covariates):
            X = sp.hstack((np.ones((self.N, 1)), covariates))
        else:
            X = covariates
        return X

    @staticmethod
    def has_bias(covariates):
        """Checks whether at least one all constant column (bias) is present in the covariates."""
        # Might have multiple invariant columns, though! Not dropped.
        for i in range(covariates.shape[1]):
            if np.all(covariates[:, i] == covariates[0, i]) and covariates[0, i] != 0:
                return True
        return False

    @staticmethod
    def _linreg(Y, X=None, Xdagger=None):
        """Efficient multiplication with symmetric covariate orthogonal projection matrix.

        Corresponds to Lippert et al., 2014, Supplement p. 10 Proposition 9
        S*a = a - X(Xdagger*a)
        with S = (I_N - X*(X.T*X)^-1*X.T) rank: N-D, symmetric covariate orthogonal projection matrix
        Source: fastlmm.association.score.py

        :param numpy.ndarray Y: factor that is multiplied with S
        :param numpy.ndarray X: covariate matrix
        :param numpy.ndarray Xdagger: precomputed Moore-Penrose pseudo-inverse of the covariate matrix or None, gets computed if not provided
        :return: result of multiplication (regressed out covariate effects); and Moore-Penrose pseudo-inverse of the covariate matrix
        :rtype: numpy.ndarray, numpy.ndarray

        The resulting term is nothing but the ordinary least squares (OLS) regression residuals after regressing out X.
        Note, that the pseudoinverse of the covariates, X†, need only be computed once (in O(ND2)) and then can be re-used across different a.
        For multiplication of matrices by S the result is applied by treating each row or column as a vector in the multiplication.
        Lippert et al., 2014

        """
        if X is None:
            RxY = Y - Y.mean(0)
            return RxY, None
        else:
            if Xdagger is None:
                Xdagger = np.linalg.pinv(X)
            RxY = Y - X.dot(Xdagger.dot(Y))
            return RxY, Xdagger

    @staticmethod
    def _hat(Y, X=None, Xdagger=None):
        """Efficient multiplication with hat-matrix.

        #TODO: add documentation
        """

        if X is None:
            Yhat = Y.mean(0)
            return Yhat, None
        else:
            if Xdagger is None:
                Xdagger = np.linalg.pinv(X)
            Yhat = X.dot(Xdagger.dot(Y))
            return Yhat, Xdagger


    def pv_alt_model(self, G1, G2=None):
        """Computes p-value of the alternative model.

        :param numpy.ndarray G1: set of SNVs to test, dimensions :math:`nxk` (:math:`n:=` number of individuals, :math:`k:=` number of SNVs to test association for)
        :return: p-value of the alternative model
        :rtype: float
        """
        if G1.shape[0] != self.N:
            logging.ERROR('Number of individuals in phenotype and genotypes to be tested does not match.')

        if G2 is None:
            squaredform, GPG = self._score(G1)
        else:
            if G1.shape[0] != G2.shape[0]:
                logging.ERROR('Number of individuals in G1 and G2 do not match.')
            squaredform, GPG = self._score_conditional(G1, G2)

        pv = self._pv_davies(squaredform, GPG)

        return pv

    def _score(self, G1):
        """Method that returns the score-based test statistic (squaredform) and the matrix (1/2)xG.TxPthetaxG (GPG) from
        which the null distribution can be efficiently computed for a set of genotypes.

        :param numpy.ndarray G1: set of SNVs to test, dimensions :math:`nxk` (:math:`n:=` number of individuals, :math:`k:=` number of SNVs to test association for)
        :return: squaredform, GPG
        :raises: NotImplementedError - Interface.
        """
        raise NotImplementedError

    def _score_conditional(self, G1, G2):
        """Method that returns the conditional score-based test statistic (squaredform) for G1 conditioned on G2 and the
        matrix (1/2)xG.TxPthetaxG (GPG) from which the null distribution can be efficiently computed for a set of genotypes.

        :param numpy.ndarray G1: set of SNVs to test, dimensions :math:`nxk` (:math:`n:=` number of individuals, :math:`k:=` number of SNVs to test association for)
        :return: squaredform, GPG
        :raises: NotImplementedError - Interface.
        """
        raise NotImplementedError

    @staticmethod
    def _pv_davies(squaredform, GPG):
        """Given the test statistic and GPG computes the corresponding p-value."""
        eigvals = LA.eigh(GPG, eigvals_only=True)
        pv = Scoretest._pv_davies_eig(squaredform, eigvals)
        return pv

    @staticmethod
    def _pv_davies_eig(squaredform, eigvals):
        """Given the test statistic and the eigenvalues of GPG computes the corresponding p-value."""
        result = Scoretest._qf(squaredform, eigvals)
        # removed keyword argument that corresponds to default.
        return result[0]

    @staticmethod
    def _qf(chi2val, coeffs, dof=None, noncentrality=None, sigma=0.0, lim=1000000, acc=1e-07):
        """Given the test statistic (squaredform) and the eigenvalues of GPG computes the corresponding p-value, calls a C script."""
        # Pia: changed default for acc to 1e-07 as this is the only value the function is ever called with in fastlmm
        from seak.cppextension import wrap_qfc

        size = coeffs.shape[0]
        if dof is None:
            dof = np.ones(size, dtype='int32')
        if noncentrality is None:
            noncentrality = np.zeros(size)
        ifault = np.zeros(1, dtype='int32')
        trace = np.zeros(7)
        pval = 1.0 - wrap_qfc.qf(coeffs, noncentrality, dof, sigma, chi2val, lim, acc, trace, ifault)
        return pval, ifault[0], trace


class ScoretestNoK(Scoretest):
    """Single kernel score-based set-association test for continuous phenotypes.

    Sets up null model for given phenotypes and covariates.
    If no covariates are given or if the covariates do not contain a bias column, a bias column is added.
    Compute p-value for alternative model with :func:`pv_alt_model`.

    Null model single kernel:

    .. math:: Y = {\\alpha} X + {\epsilon}

    Alternative model single kernel:

    .. math:: Y = {\\alpha} X + {\\gamma} G_1 + {\\epsilon}

    With: :math:`X`: covariates, dimension :math:`nxc` (:math:`n:=` number of individuals, :math:`c:=` number of covariates)
    :math:`G_1`: set of SNVs to test, dimensions :math:`nxk` (:math:`n:=` number of individuals, :math:`k:=` number of SNVs in set to test association for)

    :param numpy.ndarray phenotypes: :math:`nx1` matrix containing the phenotype (:math:`n:=` number of individuals)
    :param numpy.ndarray covariates: :math:`nxc` matrix containing the covariates (:math:`n:=` number of individuals, :math:`c:=` number of covariates)

    :ivar RxY: OLS residuals of the phenotype after regressing out fixed effects (covariates X)
    :ivar Xdagger: Moore-Penrose pseudo-inverse of the covariate matrix X
    :ivar sigma2: environmental variance
    """

    __slots__ = ["RxY", "Xdagger", "sigma2"]

    def __init__(self, phenotypes, covariates):
        """Constructor."""
        super().__init__(phenotypes, covariates)
        self.RxY, self.Xdagger, self.sigma2 = self._compute_null_model()

    def _compute_null_model(self):
        """Computes parameters of null model."""
        # residual of y regressed on X, which here, is equivalent to sigma2*Py (P is the projection matrix, which is idempotent)
        # note: Xdagger is pseudo inverse of X

        RxY, Xdagger = super()._linreg(Y=self.Y, X=self.X, Xdagger=None)

        # estimate for residual (environmental) variance
        sigma2 = (RxY * RxY).sum() / (self.Neff * self.P)

        return RxY, Xdagger, sigma2

    def _score(self, G1):
        """Computes squaredform and GPG, input for p-value computation. """

        # for the 1K case, P reduces to 1/sigma2*S
        # multiplication with S is achieved by getting the residuals regressed on X.
        # SG, needed for "GPG":
        RxG, self.Xdagger = super()._linreg(Y=G1, X=self.X, Xdagger=self.Xdagger)

        # needed for the squared form:
        GtRxY = G1.T.dot(self.RxY)

        ## original note: P is never computed explicitly, only via residuals such as Py=1/sigma2(I-Xdagger*X)y and
        ## PG=1/sigma2(I-Xdagger*X)G
        ## also note that "RxY"=Py=1/sigma2*(I-Xdagger*X)y is nothing more (except for 1/sigma2) than the residual of y
        ## regressed on X (i.e. y-X*beta), and similarly for PG="RxG"

        # note: because GtRxY has shape (D, 1), the code below is the same as (GtRxY.transpose()).dot(GtRxY)/(2 * sigma2^2):
        ## original note: yPKPy=yPG^T*GPy=(yPG^T)*(yPG^T)^T

        squaredform = ((GtRxY * GtRxY).sum()) * (0.5 / (self.sigma2 * self.sigma2))

        # we are only interested in the eigenvalues of GPG
        # np.dot(RxG.T, RxG) and np.dot(RxG, RxG.T) have the same non-zero eigenvalues!
        if G1.shape[0] > G1.shape[1]:
            # full rank, i.e. D > N
            GPG = np.dot(RxG.T, RxG)  # GPG is always a square matrix in the smaller dimension
        else:
            # low rank, i.e. D < N
            GPG = np.dot(RxG, RxG.T)

        GPG /= self.sigma2 * 2.0  # what we will take eigenvalues of for Davies, scale because P is 0.5 * 1/sigmae2 * S

        return squaredform, GPG

    def _score_conditional(self, G1, G2):
        """Computes squaredform and GPG, input for p-value computation.

        G2 are additional (fixed) effects to be conditioned on. e.g. the protein LOF burden or common variants.

        """

        # for the 1K case, P reduces to 1/sigma2*S
        # SG, needed for "GPG", i.e. GSG :
        # multiplication with S is achieved by getting the residuals regressed on X.
        # Residuals of G regressed on X:
        RxG, self.Xdagger = super()._linreg(Y=G1, X=self.X, Xdagger=self.Xdagger)

        # Residuals of G2 regressed on X:
        RxG2, self.Xdagger = super()._linreg(Y=G2, X=self.X, Xdagger=self.Xdagger)

        # we use the blockwise formula for the hat matrix, which shows that the residuals compose into 2 terms
        # when adding additional fixed effects:
        RxYhat, RxG2dagger = super()._hat(Y=self.RxY, X=RxG2, Xdagger=None)
        RxY = self.RxY + RxYhat

        # re-estimate environmental variance
        Neff = self.N - self.D - G2.shape[1]
        sigma2_update = (RxY * RxY).sum() / (Neff * self.P)

        # needed for the squared form:
        GtRxY = G1.T.dot(self.RxY) - G1.T.dot(RxG2.dot(RxG2dagger.dot(self.RxY)))

        ## original note: P is never computed explicitly, only via residuals such as Py=1/sigma2(I-Xdagger*X)y and
        ## PG=1/sigma2(I-Xdagger*X)G
        ## also note that "RxY"=Py=1/sigma2*(I-Xdagger*X)y is nothing more (except for 1/sigma2) than the residual of y
        ## regressed on X (i.e. y-X*beta), and similarly for PG="RxG"

        # note: because GtRxY has shape (D, 1), the code below is the same as (GtRxY.transpose()).dot(GtRxY)/(2 * sigma2^2):
        ## original note: yPKPy=yPG^T*GPy=(yPG^T)*(yPG^T)^T

        squaredform = ((GtRxY * GtRxY).sum()) * (0.5 / (sigma2_update * sigma2_update))

        # we are only interested in the eigenvalues of GPG
        # np.dot(RxG.T, RxG) and np.dot(RxG, RxG.T) have the same non-zero eigenvalues!
        if G1.shape[0] > G1.shape[1]:
            # full rank, i.e. D > N
            G1hat, RxG2dagger = super()._hat(G1, X=RxG2, Xdagger=RxG2dagger)
            RxG -= G1hat
            GPG = np.dot(RxG.T, RxG)  # GPG is always a square matrix in the smaller dimension
        else:
            # low rank, i.e. D < N
            G1hat, RxG2dagger = super()._hat(G1, X=RxG2, Xdagger=RxG2dagger)
            RxG -= G1hat
            GPG = np.dot(RxG, RxG.T)

        GPG /= sigma2_update * 2.0  # what we will take eigenvalues of for Davies, scale because P is 0.5 * 1/sigmae2 * S

        return squaredform, GPG


class ScoretestLogit(Scoretest):
    """Single kernel score-based set-association test for binary phenotypes.

    Sets up null model for given phenotypes and covariates.
    If no covariates are given or if the covariates do not contain a bias column, a bias column is added.
    Compute p-value for alternative model with :func:`pv_alt_model`.

    :param numpy.ndarray phenotypes: :math:`nx1` matrix containing the phenotype (:math:`n:=` number of individuals)
    :param numpy.ndarray covariates: :math:`nxc` matrix containing the covariates (:math:`n:=` number of individuals, :math:`c:=` number of covariates)
    """

    __slots__ = ["pY", "stdY", "VX", "pinvVX"]

    def __init__(self, phenotypes, covariates):
        super().__init__(phenotypes, covariates)
        # check if is binary
        uniquey = sp.unique(self.Y)
        if not sp.sort(uniquey).tolist() == [0, 1]:
            raise Exception("must use binary data in {0,1} for logit tests, found:" + str(self.Y))
        self.pY, self.stdY, self.VX, self.pinvVX = self._compute_null_model()

    def _compute_null_model(self):
        """Computes parameters of null model."""
        logreg_mod = sm.Logit(self.Y[:, 0], self.X)
        logreg_result = logreg_mod.fit(disp=0)
        pY = logreg_result.predict(self.X)
        stdY = sp.sqrt(pY * (1.0 - pY))
        VX = self.X * np.lib.stride_tricks.as_strided(stdY, (stdY.size, self.X.shape[1]), (stdY.itemsize, 0))
        pinvVX = np.linalg.pinv(VX)
        return pY, stdY, VX, pinvVX

    def _score(self, G1):
        """Computes squaredform and GPG, input for p-value computation."""
        RxY = (self.Y.flatten() - self.pY)  # residual of y regressed on X, which here, is equivalent to sigma2*Py
        # (P is the projection matrix, which is idempotent)
        VG = G1 * np.lib.stride_tricks.as_strided(self.stdY, (self.stdY.size, G1.shape[1]), (self.stdY.itemsize, 0))
        GY = G1.T.dot(RxY)
        squaredform = (GY * GY).sum() / (2.0 * self.P)
        RxVG, Xd = super()._linreg(VG, X=self.VX, Xdagger=self.pinvVX)

        if G1.shape[0] < G1.shape[1]:
            GPG = RxVG.dot(RxVG.T) / (2.0 * self.P)
        else:
            GPG = RxVG.T.dot(RxVG) / (2.0 * self.P)
        return squaredform, GPG


class Scoretest2K(Scoretest):
    """Two-kernel score-based set-association test for continuous phenotypes.

    Sets up null model for given phenotypes, covariates and background kernel :math:`K_0` or background genotypes :math:`G_0`.
    If no covariates are given or if the covariates do not contain a bias column, a bias column is added.
    Compute p-value for alternative model with :func:`pv_alt_model`.

    Null model two-kernel:

    .. math:: Y = {\\alpha} X + {\\beta} G_0 + {\\epsilon}

    Alternative model two-kernel:

    .. math:: Y = {\\alpha} X + {\\beta} G_0 + {\\gamma} G_1 + {\\epsilon}

    With: :math:`X`: covariates, dimension :math:`nxc` (:math:`n:=` number of individuals, :math:`c:=` number of covariates)
    :math:`G_0`: SNVs to construct the background kernel :math:`K_0` from, correcting for population structure, dimensions :math:`nxm` (:math:`n:=` number of individuals, :math:`m:=` number of SNVs/variants)
    :math:`G_1`: set of SNVs to test, dimensions :math:`nxk` (:math:`n:=` number of individuals, :math:`k:=` number of SNVs in set to test association for)

    :param numpy.ndarray phenotypes: :math:`nx1` matrix containing the phenotype (:math:`n:=` number of individuals)
    :param numpy.ndarray covariates: :math:`nxc` matrix containing the covariates (:math:`n:=` number of individuals, :math:`c:=` number of covariates)
    :param numpy.ndarray K0: genetic similarity matrix/GRM :math:`K_0` accounting for confounding
    :param numpy.ndarray G0: genotype matrix :math:`G_0` used to contruct math:`K_0` to account for confounding
    :param boolean forcefullrank: for testing purposes only

    :ivar S: eigenvalues of PKP
    :ivar Xdagger: Moore-Penrose pseudo-inverse of the covariate matrix X
    :ivar sigma2e: environmental variance
    :ivar sigma2g: genetic variance
    """

    __slots__ = ["K0", "G0", "Xdagger", "S", "U", "lowrank", "UY", "UUY", "YUUY", "optparams", "sigma2e", "sigma2g"]

    def __init__(self, phenotypes, covariates=None, K0=None, G0=None, forcefullrank=False):
        """Constructor."""
        # note: super().__init__ simply fills the slots for covariates etc. no computation done yet:
        super().__init__(phenotypes, covariates)
        self.G0 = G0
        self.K0 = K0
        # spectral decomposition, needed to efficiently compute the matrix square root of P, see Lippert 2014, suppl. 7.3
        self.Xdagger, self.S, self.U, self.lowrank, self.UY, self.YUUY = self._compute_spectral_decomposition(forcefullrank=forcefullrank)
        self.optparams = self._compute_null_model()
        self.sigma2e = (1.0 - self.optparams["h2"]) * self.optparams["sigma2"]
        self.sigma2g = self.optparams["h2"] * self.optparams["sigma2"]

    def _compute_spectral_decomposition(self, forcefullrank=False):
        """Computes spectral decomposition of K0 or G0."""
        lowrank = False
        Xdagger = None

        # we want to compute the spectral decomposition of P (needed for matrix square root)
        # sigma2g * P = (S(Kg+delta*I)S)^dagger
        # we can compute the spectral decomposition of the right hand


        if self.K0 is not None:
            # K0 already computed, i.e. the full rank case, work with K0
            # compute SVD of SKS
            # see Lemma 10 in Lippert et al. 2014
            ar = sp.arange(self.K0.shape[0])

            self.K0[ar, ar] += 1.0

            PxKPx, Xdagger = super()._linreg(Y=self.K0, X=self.X, Xdagger=Xdagger)
            PxKPx, self.Xdagger = super()._linreg(Y=PxKPx.T, X=self.X, Xdagger=Xdagger)

            # S are the eigenvalues
            # U are the eigenvectors
            [S, U] = LA.eigh(PxKPx)
            self.K0[ar, ar] -= 1.0
            U = U[:, self.D:self.N]
            S = S[self.D:self.N] - 1.0

        elif 0.7 * (self.Neff) <= self.G0.shape[1] or forcefullrank:
            # K0 gets computed, work with K0

            self.K0 = self.G0.dot(self.G0.T)
            # compute SVD of SKS
            # see Lemma 10 in Lippert et al. 2014
            # the rest is identical to the case above...
            ar = sp.arange(self.K0.shape[0])
            self.K0[ar, ar] += 1.0

            PxKPx, Xdagger = super()._linreg(Y=self.K0, X=self.X, Xdagger=Xdagger)
            PxKPx, self.Xdagger = super()._linreg(Y=PxKPx.T, X=self.X, Xdagger=Xdagger)

            # S are the eigenvalues
            # U are the eigenvectors
            self.K0[ar, ar] -= 1.0
            [S, U] = LA.eigh(PxKPx)
            U = U[:, self.D:self.N]
            S = S[self.D:self.N] - 1.0

        else:
            # work with G0 instead of K0, this is the low-rank case

            PxG, Xdagger = super()._linreg(Y=self.G0, X=self.X, Xdagger=Xdagger)
            [U, S, V] = LA.svd(PxG, False, True)
            inonzero = S > 1E-10
            # S are the eigenvalues
            # U are the eigenvectors
            S = S[inonzero] * S[inonzero]
            U = U[:, inonzero]
            lowrank = True

        UY = U.T.dot(self.Y)

        if lowrank:
            Yres, Xdagger = super()._linreg(Y=self.Y, X=self.X, Xdagger=Xdagger)
            UUY = Yres - U.dot(UY)
            YUUY = (UUY * UUY).sum()
        else:
            YUUY = None

        return Xdagger, S, U, lowrank, UY, YUUY

    def _compute_null_model(self):
        """Computes parameters of null model."""
        resmin = [None]

        def f(x, resmin=resmin, **kwargs):
            res = self._nLLeval(h2=x)
            if (resmin[0] is None) or (res['nLL'] < resmin[0]['nLL']):
                resmin[0] = res
            return res['nLL']

        minimize1D(f, evalgrid=None, nGrid=20, minval=0.0, maxval=0.99999)

        # dictionary containing the model parameters at the optimal h2
        optparams = resmin[0]

        return optparams

    def _nLLeval(self, h2=0.0):
        """
        evaluate -ln( N( U^T*y | U^T*X*beta , h2*S + (1-h2)*I ) ),
        where K = USU^T
        --------------------------------------------------------------------------
        Input:
        h2      : mixture weight between K and Identity (environmental noise)
        --------------------------------------------------------------------------
        Output dictionary:
        'nLL'       : negative log-likelihood
        'sigma2'    : the model variance sigma^2
        'h2'        : mixture weight between Covariance and noise
        --------------------------------------------------------------------------
        """
        if (h2 < 0.0) or (h2 >= 1.0):
            return {'nLL': 3E20,
                    'h2': h2
                    }
        k = self.S.shape[0]

        Sd = h2 * self.S + (1.0 - h2)
        UYS = self.UY / np.lib.stride_tricks.as_strided(Sd, (Sd.size, self.UY.shape[1]), (Sd.itemsize, 0))

        YKY = (UYS * self.UY).sum()

        logdetK = sp.log(Sd).sum()

        if (self.lowrank):  # low rank part
            YKY += self.YUUY / (1.0 - h2)
            logdetK += sp.log(1.0 - h2) * (self.Neff * self.P - k)

        sigma2 = YKY / (self.Neff * self.P)
        nLL = 0.5 * (logdetK + self.Neff * self.P * (sp.log(2.0 * sp.pi * sigma2) + 1))
        result = {
            'nLL': nLL,
            'sigma2': sigma2,
            'h2': h2
        }
        return result


    def _score(self, G1):
        """Computes squaredform and GPG with a background kernel."""

        # SG
        RxG, Xdagger = super()._linreg(Y=G1, X=self.X, Xdagger=self.Xdagger)

        # UtSG
        UG = self.U.T.dot(RxG)

        if self.lowrank:
            UUG = RxG - self.U.dot(UG)

        # Compare to Lippert 2014, Lemma 11. Rescale eigenvalues according to mixing parameters
        # The inverse of the diagonal matrix of eigenvalues (Lambda + delta*I) is calculated trivially:
        Sd = 1.0 / (self.S * self.sigma2g + self.sigma2e)

        # matrix multiplication of UtSG with (Lambda + delta*I)^-1, which is called Sd here
        SUG = UG * np.lib.stride_tricks.as_strided(Sd, (Sd.size, UG.shape[1]), (Sd.itemsize, 0))

        GPY = SUG.T.dot(self.UY)
        if self.lowrank:
            GPY += UUG.T.dot(self.UUY) / self.sigma2e

        # see GPY in (14)
        squaredform = 0.5 * (GPY * GPY).sum()

        if G1.shape[0] > G1.shape[1]:
            GPG = SUG.T.dot(UG)
        else:
            GPG = SUG.dot(UG.T)

        # in the original they compute expectationsqform (expected value) and varsqform (variance) fo the squared form.
        # these were not used.

        if self.lowrank:
            if G1.shape[0] > G1.shape[1]:
                GPG_lowr = UUG.T.dot(UUG) / self.sigma2e
            else:
                GPG_lowr = UUG.dot(UUG.T) / self.sigma2e
            GPG += GPG_lowr

        GPG = GPG * 0.5

        return squaredform, GPG
