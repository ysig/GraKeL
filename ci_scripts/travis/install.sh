# Configure the conda environment and put it in the path using the
# provided versions
if [[ "$TRAVIS_OS_NAME" == "osx" ]]; then
    brew update;
    brew install openssl;
    ln -s /usr/local/opt/openssl/lib/libcrypto.1.0.0.dylib /usr/local/lib/
    ln -s /usr/local/opt/openssl/lib/libssl.1.0.0.dylib /usr/local/lib/
    ln -s /usr/local/Cellar/openssl/1.0.2j/bin/openssl /usr/local/bin/openssl
    brew reinstall {python} --with-brewed-openssl
    brew link --overwrite {python}
fi

{pip} install --upgrade pip
{pip} install --upgrade setuptools
{pip} install -r requirements.txt
{pip} install "cvxopt==1.2.0"

# Install coverage if needed
if [[ "$COVERAGE" == "true" ]]; then
    {pip} install coverage
fi

# Print versions
{python} --version
{python} -c "import numpy; print('numpy %s' % numpy.__version__)"
{python} -c "import scipy; print('scipy %s' % scipy.__version__)"
{python} -c "import cython; print('cython %s' % cython.__version__)"
{python} -c "import cvxopt; print('cvxopt %s' % cvxopt.__version__)"

