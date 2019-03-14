# !/bin/bash

make clean

make -j8

make codecs

rm test1.log

#./codecs --codecs pfor --needtodelta --zipfian1
./codecs --codecs iiu --needtodelta --needtoskiplist --zipfian1 
#./codecs --codecs pfor,iiu --needtodelta --zipfian1
