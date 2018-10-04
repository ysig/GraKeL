set -e
current_dir="$(pwd)"

if [[ "$COVERAGE" == "true" ]]; then
    # Store artifacts
    mkdir -p $TEST_DIR;
    cd $TEST_DIR;
    pwd;
    nosetests --with-coverage --cover-package=$MODULE $MODULE;
    ls -x -l;
else
    # Ignore arifacts: just change folder
    mkdir -p $TEST_DIR;
    cd $TEST_DIR;
    nosetests $MODULE;
fi

cd $current_dir
