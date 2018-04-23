#!/bin/bash
# This script is meant to be called in the "deploy" step defined in 
# circle.yml. See https://circleci.com/docs/ for more details.
# The behavior of the script is controlled by environment variable defined
# in the circle.yml in the top level folder of the project.

MSG="Pushing the docs for revision for branch: $CIRCLE_BRANCH, commit $CIRCLE_SHA1"

cd $HOME
# Copy the build docs to a temporary folder

# rename the project folder to the doc repo folder
if [ -d project ];
    then mv project $DOC_REPO;
fi

rm -rf tmp
mkdir tmp
cp -R $HOME/$DOC_REPO/doc/_build/html/* ./tmp/ 

# Clone the docs repo if it isnt already there
if [ ! -d $DOC_REPO ];
    then git clone "git@github.com:$USERNAME/"$DOC_REPO".git";
fi

cd $DOC_REPO
git branch gh-pages
git checkout -f gh-pages
git reset --hard origin/gh-pages
git clean -dfx

for name in $(ls -A $HOME/$DOC_REPO); do
    case $name in
        .nojekyll) # So that github does not build this as a Jekyll website.
        ;;
        circle.yml) # Config so that build gh-pages branch.
        ;;
        *)
        git rm -rf $name
        ;;
    esac
done

# Copy the new build docs
mkdir $DOC_URL
cp -R $HOME/tmp/* ./$DOC_URL/

git config --global user.email $EMAIL
git config --global user.name $USERNAME
git add -f ./$DOC_URL/
git commit -m "$MSG"
git push -f origin gh-pages
if [ $? -ne 0 ]; then
    echo "Pushing docs failed"
    echo
    exit 1
fi

echo $MSG

# Deploy on PyPi
if [[ $DEPLOY_PYPI == "true" ]]; then
    # Move to home repo
    cd $HOME/$DOC_REPO
    ls -x

    # Initialise .pypirc
    echo "[distutils]" > ~/.pypirc
    echo "index-servers = pypi" >> ~/.pypirc
    echo "[pypi]" >> ~/.pypirc
    echo "repository=https://pypi.python.org/pypi" >> ~/.pypirc
    echo "username=$USERNAME" >> ~/.pypirc
    echo "password=$PYPI_PASSWORD" >> ~/.pypirc

    # Initialise setup.cfg
    echo "[build_sphinx]" >> .setup.cfg
    echo "source-dir = doc" >> .setup.cfg
    echo "build-dir  = doc/_build" >> .setup.cfg
    echo "all_files  = 1" >> .setup.cfg
    echo "[upload_sphinx]" >> .setup.cfg
    echo "upload-dir = doc/_build/html" >> .setup.cfg

    # Upload source files
    rm -rf dist/
    sudo pip install twine --upgrade
    sudo python setup.py sdist --formats=zip
    twine upload dist/*.zip

    # Build & Upload sphinx-docs
    sudo pip install sphinx-pypi-upload
    sudo python setup.py build_sphinx
    sudo python setup.py upload_sphinx
fi
