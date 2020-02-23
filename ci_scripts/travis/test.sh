set -e
current_dir="$(pwd)"

if [[ "$COVERAGE" == "true" ]]; then
    # Store artifacts
    pip install coverage
    pip install codecov
    nosetests $MODULE --with-coverage --cover-package=$MODULE;

    # Ignore coveralls failures as the coveralls server is not
    # very reliable but we don't want travis to report a failure
    # in the github UI just because the coverage report failed to
    # be published.
    cp -r .coverage /project/.coverage
    cd /project
    echo "Coverage";
    ls -la;
    python -m codecov;
else
    # Ignore arifacts: just change folder
    mkdir -p $TEST_DIR;
    cd $TEST_DIR;
    nosetests $MODULE;
fi

cd $current_dir
