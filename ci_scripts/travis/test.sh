set -e
current_dir="$(pwd)"

if [[ "$COVERAGE" == "true" ]]; then
    # Store artifacts
    mkdir -p $TEST_DIR;
    cd $TEST_DIR;
    pwd;
    nosetests $MODULE --with-coverage --cover-package=$MODULE --cover-xml;
    
    echo "Coverage"
    ls -la;
  
    # Ignore coveralls failures as the coveralls server is not
    # very reliable but we don't want travis to report a failure
    # in the github UI just because the coverage report failed to
    # be published.
    $PYTHON -m codecov
fi
else
    # Ignore arifacts: just change folder
    mkdir -p $TEST_DIR;
    cd $TEST_DIR;
    nosetests $MODULE;
fi

cd $current_dir
