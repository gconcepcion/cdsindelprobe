# -*- coding: utf-8 -*-
import os
from setuptools import setup, find_packages


_REQUIREMENTS_FILE = 'REQUIREMENTS.txt'

with open('README.rst') as f:
    readme = f.read()

with open('LICENSE') as f:
    license = f.read()

def _get_local_file(file_name):
    return os.path.join(os.path.dirname(__file__), file_name)

def _get_requirements(file_name):
    with open(file_name, 'r') as f:
        lines = f.readlines()
    reqs = [l for l in lines if not l.startswith("#")]
    return reqs

setup(
    name='cdsindelprobe',
    version='0.0.1',
    description='Summarize 1 bp frameshifts in CDS aligned to reference genome',
    long_description=readme,
    author='Greg Concepcion',
    author_email='gconcepcion@pacificbiosciences.com',
    url='https://github.com/gconcepcion/cdsindelprobe',
    license=license,
    scripts=['bin/probe_indels.py'],
    packages=find_packages(),
    zip_safe=False,
    install_requires=_get_requirements(_get_local_file(_REQUIREMENTS_FILE))
)

