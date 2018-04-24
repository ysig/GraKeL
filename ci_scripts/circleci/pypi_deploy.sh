#!/bin/bash
# This script is meant to be called in the "deploy" step defined in 
# circle.yml. See https://circleci.com/docs/ for more details.
# The behavior of the script is controlled by environment variable defined
# in the circle.yml in the top level folder of the project.

# Deploy on PyPi
if [[ $DEPLOY_PYPI == "true" ]]; then
    # Build & Upload sphinx-docs
    . venv/bin/activate
    python setup.py upload_sphinx || true
fi
