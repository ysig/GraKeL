set -e
current_dir="$(pwd)"

if [[ "$COVERAGE" == "true" ]]; then
    # Store artifacts
    ls /host
    cp -r $TRAVIS_BUILD_DIR/.git .
    cp  $TRAVIS_BUILD_DIR/git .

    pip install coverage
    pip install codecov
    nosetests $MODULE --with-coverage --cover-package=$MODULE;

    echo "Coverage";
    ls -la;

    # Ignore coveralls failures as the coveralls server is not
    # very reliable but we don't want travis to report a failure
    # in the github UI just because the coverage report failed to
    # be published.
    python -m codecov;
else
    # Ignore arifacts: just change folder
    mkdir -p $TEST_DIR;
    cd $TEST_DIR;
    nosetests $MODULE;
fi

cd $current_dir
