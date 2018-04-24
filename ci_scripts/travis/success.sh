set -e
    
if [[ "$COVERAGE" == "true" ]] && [[ "$TRAVIS_OS_NAME" == "linux" ]] && [[ "$TRAVIS_PYTHON_VERSION" == "2.7" ]]; then
    # Need to run coveralls from a git checkout, so we copy .coverage
    # from TEST_DIR where nosetests has been run
    cp $TEST_DIR/.coverage $TRAVIS_BUILD_DIR
    cd $TRAVIS_BUILD_DIR
    # Ignore coveralls failures as the coveralls server is not
    # very reliable but we don't want travis to report a failure
    # in the github UI just because the coverage report failed to
    # be published.
    coveralls || echo "Coveralls upload failed"
fi

if [[ "$DEPLOY_WHEEL" == "true" ]]; then
    pip show cython
    pip install twine --upgrade
fi

if [[ "$DEPLOY_WHEEL" == "true" ]] && [[ "$TRAVIS_OS_NAME" == "linux" ]] && [[ "$TRAVIS_PYTHON_VERSION" == "2.7" ]]; then
    # Upload source files
    rm -rf dist/
    python setup.py sdist --formats=zip
    twine upload dist/*.zip || true
fi

if [[ "$DEPLOY_WHEEL" == "true" ]]; then
    cd $TRAVIS_BUILD_DIR
    pip install cibuildwheel==0.7.1
    cibuildwheel --output-dir wheelhouse
    twine upload wheelhouse/*.whl || true
fi
