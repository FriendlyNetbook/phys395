#!/bin/bash

gfortran -O3 fibonacci.f90 

./a.out > DATA

gnuplot -e "plot 'DATA'" 

