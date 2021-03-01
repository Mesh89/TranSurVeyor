tar -xjf htslib-1.11.tar.bz2
cd htslib-1.11
autoheader
autoconf
./configure --prefix=`pwd`
make
make install
