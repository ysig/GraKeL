#!/bin/bash
# This script is meant to be called in the "deploy" step defined in 
# circle.yml. See https://circleci.com/docs/ for more details.
# The behavior of the script is controlled by environment variable defined
# in the circle.yml in the top level folder of the project.

# System dependencies
sudo -E apt-get -yq remove texlive-binaries --purge > /dev/null
sudo apt-get update > /dev/null
sudo apt-get install libatlas-dev libatlas3gf-base > /dev/null
sudo apt-get install build-essential python-dev python-setuptools > /dev/null

# Setup a python venv and install basics
source ~/venv/bin/activate
pip install --upgrade pip
pip install --upgrade matplotlib setuptools nose coverage sphinx pillow sphinx-gallery sphinx_rtd_theme sphinxcontrib-bibtex nb2plots numpydoc tqdm > /dev/null
pip install -r requirements.txt > /dev/null
pip install "cvxopt==1.2.0" > /dev/null

# More dependencies
sudo -E apt-get -yq update > /dev/null
sudo -E apt-get -yq --no-install-suggests --no-install-recommends --force-yes install dvipng texlive-latex-base texlive-latex-extra > /dev/null

# Install project
python setup.py clean
python setup.py develop

# Build Docs
set -o pipefail && cd doc && make html 2>&1 | tee ~/log.txt
cat ~/log.txt && if grep -q "Traceback (most recent call last):" ~/log.txt; then false; else true; fi

# Build for deploy
if [[ $DEPLOY_PYPI == "true" ]]; then
    # Initialise .pypirc
    echo "[distutils]" > ~/.pypirc
    echo "index-servers = pypi" >> ~/.pypirc
    echo "[pypi]" >> ~/.pypirc
    echo "username=$USERNAME" >> ~/.pypirc
    echo "password=$PYPI_PASSWORD" >> ~/.pypirc

    # Initialise setup.cfg
    echo "[build_sphinx]" >> .setup.cfg
    echo "source-dir = doc" >> .setup.cfg
    echo "build-dir  = doc/_build" >> .setup.cfg
    echo "all_files  = 1" >> .setup.cfg
    echo "[upload_sphinx]" >> .setup.cfg
    echo "upload-dir = doc/_build/html" >> .setup.cfg

    # Build sphinx docs
    pip install sphinx-pypi-upload
    python setup.py build_sphinx
fi
