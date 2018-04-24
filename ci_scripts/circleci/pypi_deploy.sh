#!/bin/bash
# This script is meant to be called in the "deploy" step defined in 
# circle.yml. See https://circleci.com/docs/ for more details.
# The behavior of the script is controlled by environment variable defined
# in the circle.yml in the top level folder of the project.

# Deploy on PyPi
if [[ $DEPLOY_PYPI == "true" ]]; then
    # Download clean repo
    rm -rf ~/tmp && mkdir ~/tmp && cd ~/tmp
    git clone "git@github.com:$USERNAME/"$DOC_REPO".git"
    cd $DOC_REPO
    

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

    # Build & Upload sphinx-docs
    . venv/bin/activate
    pip install sphinx-pypi-upload
    python setup.py build_sphinx
    python setup.py upload_sphinx || true
fi
