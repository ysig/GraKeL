# Install dependencies related with the [SSL: TLSV1_ALERT_PROTOCOL_VERSION]
# a urllib2.URLError
if [[ "$TRAVIS_OS_NAME" == "osx" ]]; then
    echo "normal openssl";
    openssl version -a;
    {python} -c "import ssl; print('py:openssl:', ssl.OPENSSL_VERSION)";
    brew install {python} --with-brewed-openssl
    {python} -c "import ssl; print('py:openssl:', ssl.OPENSSL_VERSION)";   
#    brew unlink {python}
#    brew reinstall {python} --with-brewed-openssl
#    brew link --overwrite {python}
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

