#!/usr/bin/env sh

# needs root privileges

# Requirements
apt-get install python-dev
apt-get -y install python3-pip
pip3 install --upgrade setuptools

# Download library
wget https://web.cs.dal.ca/~peter/software/pynauty/pynauty-0.6.0.tar.gz

# decompress
tar -xf pynauty-0.6.0.tar.gz && rm -f pynauty-0.6.0.tar.gz

# change directory
cd pynauty-0.6.0

# download c library
wget http://users.cecs.anu.edu.au/~bdm/nauty/nauty26r10.tar.gz
# decompress library
tar -xf nauty26r10.tar.gz

# make a soft link
ln -s nauty26r10 nauty

# build pynauty
make pynauty

ifdef VIRTUAL_ENV
	make virtenv-ins
else
	make tests
	make user-ins
endif

# exit directory and delete
cd ..
rm -rf pynauty-0.6.0
