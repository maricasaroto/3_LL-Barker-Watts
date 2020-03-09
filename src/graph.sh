#!/bin/bash
 
# created     : 2019/05/15
# last update : 2019/07/22
# author      : Mariana Casaroto <mariana.fcasaroto@gmail.com>
# notes       : confecção dos gráficos

cd ../graph/S/
gnuplot S.plt
make
cd ../E/
gnuplot E.plt
make
cd ../acc/
gnuplot acc.plt
make
cd ../sig/
gnuplot sig.plt
make
cd ../binder/
gnuplot binder.plt
make
cd ../E2/
gnuplot E2.plt
make
cd ../E4/
gnuplot E4.plt
make
cd ../sigS/
gnuplot sigS.plt
make
cd ../ref/
gnuplot ref.plt
make
cd ../Cv/
gnuplot cv.plt
make

