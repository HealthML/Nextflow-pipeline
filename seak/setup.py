#!/usr/bin/env python
# -*- encoding: utf-8 -*-
from __future__ import absolute_import
from __future__ import print_function

import io
import re
from glob import glob
from os.path import basename
from os.path import dirname
from os.path import join
from os.path import splitext

import numpy
from setuptools import Extension
from setuptools import find_packages
from setuptools import setup


def read(*names, **kwargs):
    with io.open(
        join(dirname(__file__), *names),
        encoding=kwargs.get('encoding', 'utf8')
    ) as fh:
        return fh.read()


ext_modules = [Extension(name="seak.cppextension.wrap_qfc",
                         language="c++",
                         sources=["src/seak/cppextension/wrap_qfc.cpp", "src/seak/cppextension/QFC.cpp"],
                         include_dirs=[numpy.get_include()])]  # , define_macros=macros)]

setup(
    name='seak',
    version='0.0.0',
    license='Apache-2.0',
    description='Program that performs kernel-based genome-wide association studies of sets of SNVs.',
    long_description='%s\n%s' % (
        re.compile('^.. start-badges.*^.. end-badges', re.M | re.S).sub('', read('README.rst')),
        re.sub(':[a-z]+:`~?(.*?)`', r'``\1``', read('CHANGELOG.rst'))
    ),
    author='Pia Rautenstrauch',
    author_email='None',
    url='https://github.com/rautensp/seak_thesis_submission',
    packages=find_packages('src'),
    package_dir={'': 'src'},
    py_modules=[splitext(basename(path))[0] for path in glob('src/*.py')],
    include_package_data=True,
    zip_safe=False,
    classifiers=[
        # complete classifier list: http://pypi.python.org/pypi?%3Aaction=list_classifiers
        'Development Status :: 5 - Production/Stable',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: Apache Software License',
        'Operating System :: Unix',
        'Operating System :: POSIX',
        'Operating System :: Microsoft :: Windows',
        'Programming Language :: Python',
        # 'Programming Language :: Python :: 2.7',
        # 'Programming Language :: Python :: 3',
        # 'Programming Language :: Python :: 3.5',
        # 'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        # 'Programming Language :: Python :: Implementation :: CPython',
        # 'Programming Language :: Python :: Implementation :: PyPy',
        # uncomment if you test on these interpreters:
        # 'Programming Language :: Python :: Implementation :: IronPython',
        # 'Programming Language :: Python :: Implementation :: Jython',
        # 'Programming Language :: Python :: Implementation :: Stackless',
        'Topic :: Utilities',
        'Private :: Do Not Upload',
    ],
    project_urls={
        # 'Documentation': 'https://seak.readthedocs.io/',
        # 'Changelog': 'https://seak.readthedocs.io/en/latest/changelog.html',
        'Issue Tracker': 'https://github.com/rautensp/seak_thesis_submission/issues',
    },
    keywords=[
        # eg: 'keyword1', 'keyword2', 'keyword3',
    ],
    python_requires='>=2.7, !=3.0.*, !=3.1.*, !=3.2.*, !=3.3.*, !=3.4.*',
    install_requires=['numpy>=1.17.2', 'pandas>=0.25.3', 'h5py', 'pandas_plink', 'pyreadr', 'scipy', 'statsmodels',
                      'matplotlib', 'scikit-learn>=0.21', 'pysnptools'  # eg: 'aspectlib==1.1.1', 'six>=1.7',
                      ],
    extras_require={
        # eg:
        #   'rst': ['docutils>=0.11'],
        #   ':python_version=="2.6"': ['argparse'],
    },
    entry_points={
        'console_scripts': [
            'seak = seak.cli:main',
        ]
    },
    ext_modules=ext_modules,
)
