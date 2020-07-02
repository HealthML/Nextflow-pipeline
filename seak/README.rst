========
Overview
========

This documentation portrays the status quo of the package :mod:`seak` at the time of the submission of the Master's thesis of Pia Rautenstrauch (31.03.2020).

:mod:`Seak`, which stands for **se**\ quence **a**\ nnotations in **k**\ ernel-based tests, is an open-source Python
software package for performing set-based genotype-phenotype association tests. It allows for the flexible incorporation
of prior knowledge, such as variant effect predictions, or other annotations, into variant association tests via kernel
functions.  The mathematical implementation of these tests is based on
:mod:`FaST-LMM-Set` :cite:`Listgarten2013` :cite:`Lippert2014`. It can correct for population and family structure as well as
cryptic relatedness in a two random effects model (:class:`seak.scoretest.Scoretest2K`) for continuous phenotypes.
Furthermore, it also implements a one random effect model with a linear (:class:`seak.scoretest.ScoretestNoK`)
and a logistic (:class:`seak.scoretest.ScoretestLogit`) link function. Associations are tested in variance component score tests.

In addition to the supported input file formats, I implemented interfaces for all data loading functionalities (:mod:`seak.data_loaders`).
This way, I seek to maximize the flexibility of :mod:`seak`, such that expert users can easily adapt the package to the input data types of their choice.

* Free software: Apache Software License 2.0

Installation
============
The installation of :mod:`seak` requires Python 3.7+ and the packages `numpy <https://pypi.org/project/numpy/>`_ and `cython <https://pypi.org/project/Cython/>`_. All other dependencies are installed automatically when installing the package.

At the command line::

    python seak_thesis_submission/setup.py install


Documentation
=============
For a reference documenting all public modules included in :mod:`seak` meant for general usage see:
:ref:`API reference`.

Tutorial
========
A small example illustrating how to perform a two variance component score test with :mod:`seak` is shown in: :ref:`Tutorial`.


References
=============

.. bibliography:: references.bib
    :style: unsrt
