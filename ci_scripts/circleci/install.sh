#!/bin/bash
# This script is meant to be called in the "deploy" step defined in 
# circle.yml. See https://circleci.com/docs/ for more details.
# The behavior of the script is controlled by environment variable defined
# in the circle.yml in the top level folder of the project.

# System dependencies
sudo -E apt-get -yq remove texlive-binaries --purge
sudo apt-get update
sudo apt-get install libatlas-dev libatlas3gf-base
sudo apt-get install build-essential python-dev python-setuptools

# Setup a python venv and install basics
python -m venv venv
. venv/bin/activate
pip install --upgrade numpy
pip install --upgrade scipy matplotlib setuptools nose coverage sphinx pillow sphinx-gallery sphinx_rtd_theme sphinxcontrib-bibtex nb2plots

# More dependencies
sudo -E apt-get -yq update
sudo -E apt-get -yq --no-install-suggests --no-install-recommends --force-yes install dvipng texlive-latex-base texlive-latex-extra

# Post dependencies
pip install --upgrade tqdm
pip install --upgrade cython numpydoc
pip install --upgrade scikit-learn
pip install --upgrade "cvxopt>=1.2.0"

# Install project
python setup.py clean
python setup.py develop

# Build Docs
set -o pipefail && cd doc && make html 2>&1 | tee ~/log.txt
cat ~/log.txt && if grep -q "Traceback (most recent call last):" ~/log.txt; then false; else true; fi
