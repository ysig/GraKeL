set -e
current_dir="$(pwd)"

if [[ "$COVERAGE" == "true" ]]; then
    # Store artifacts: user is mounted on host
    if [[ "$TRAVIS_OS_NAME" == "linux" ]]; then
        mkdir -p "/host$TEST_DIR";
        cd "/host$TEST_DIR";
    else
        mkdir -p $TEST_DIR;
        cd $TEST_DIR;
    fi
    pwd;
    ls;
    nosetests --with-coverage --cover-package=$MODULE $MODULE;
else
    # Ignore arifacts: just change folder
    mkdir -p $TEST_DIR;
    cd $TEST_DIR;
    nosetests $MODULE;
fi

cd $current_dir
