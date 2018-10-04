set -e
current_dir="$(pwd)"

if [[ "$COVERAGE" == "true" ]]; then
    # Store artifacts
    mkdir -p $TEST_DIR;
    cd $TEST_DIR;
    pwd;
    nosetests $MODULE --with-coverage --cover-package=$MODULE --cover-xml;
    ls -la;
else
    # Ignore arifacts: just change folder
    mkdir -p $TEST_DIR;
    cd $TEST_DIR;
    nosetests $MODULE;
fi

cd $current_dir
