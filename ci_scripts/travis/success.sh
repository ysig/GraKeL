set -e
cd $TRAVIS_BUILD_DIR
rm -rf dist/

$PIP install twine
$PIP install setuptools --upgrade

if [[ "$DEPLOY_SDIST" == "true" ]]; then
    # Build & Deploy sdist
    $PYTHON setup.py sdist --formats=zip 
    $PYTHON -m twine upload dist/*.zip || true
fi

if [[ "$DEPLOY_WHEEL" == "true"  ]]; then
    # Deploy wheels
    $PYTHON -m twine upload $WHEEL_FOLDER/*.whl || true
fi
