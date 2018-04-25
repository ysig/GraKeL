#!/bin/bash
# This script is meant to be called in the "deploy" step defined in 
# circle.yml. See https://circleci.com/docs/ for more details.
# The behavior of the script is controlled by environment variable defined
# in the circle.yml in the top level folder of the project.

# Deploy on PyPi: Only works with python2
if [[ $DEPLOY_PYPI == "true" ]]; then
    # Build & Upload sphinx-docs
    # Initialise .pypirc
    echo "[distutils]" > ~/.pypirc
    echo "index-servers = pypi" >> ~/.pypirc
    echo >> ~/.pypirc
    echo "[pypi]" >> ~/.pypirc
    echo "username=$USERNAME" >> ~/.pypirc
    echo "password=$PYPI_PASSWORD" >> ~/.pypirc
    
    # Upload sphinx docs
    ls ~/project
    ls ~/project/doc
    cd ~/project
    sudo apt-get install tree
    tree -d ~/project || true
    source ~/project/venv/bin/activate
    pip install sphinx-pypi-upload
    python setup.py upload_sphinx --upload-dir=~/project/doc/_build/html || true
fi
