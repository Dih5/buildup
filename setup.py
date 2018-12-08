#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""The setup script."""

from setuptools import setup, find_packages

with open('README.rst') as readme_file:
    readme = readme_file.read()

with open('HISTORY.rst') as history_file:
    history = history_file.read()

requirements = ['physdata']

setup_requirements = ['pytest-runner', ]

test_requirements = ['pytest']

setup(
    author="Dih5",
    author_email='dihedralfive@gmail.com',
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: GNU Lesser General Public License v3 or later (LGPLv3+)',
        'Natural Language :: English',
        "Programming Language :: Python :: 2",
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
    ],
    description="Build up factor",
    entry_points={
        'console_scripts': [
        ],
        'gui_scripts': [
        ],

    },
    install_requires=requirements,
    license='LGPLv3+',
    long_description=readme + '\n\n' + history,
    include_package_data=True,
    keywords=['build up', 'spectrometry', 'particle transport'],
    name='buildup',
    package_data={
        'buildup': ['buildup/data','buildup/model'],
    },
    packages=find_packages(include=['buildup']),
    setup_requires=setup_requirements,
    test_suite='tests',
    tests_require=test_requirements,
    url='https://github.com/Dih5/buildup',
    version='0.1.0',
    zip_safe=False,
)
