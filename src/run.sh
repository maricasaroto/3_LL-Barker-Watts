#!/bin/bash
 
# created     : 2019/05/15
# last update : 2019/05/15
# author      : Mariana Casaroto <mariana.fcasaroto@gmail.com>
# notes       : 
 
gcc -O3 -lm -lgsl -lgslcblas bw.c 
time ./a.out
#./graph.sh
