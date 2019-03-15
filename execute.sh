# !/bin/bash

make clean

make -j8

make codecs

rm test1.log

#./codecs --codecs pfor --needtodelta --zipfian1
#./codecs --codecs pfor,iiu --needtodelta --zipfian1
#./codecs --codecs pfor,iiu --needtodelta --needtoskiplist 
#./codecs --codecs pfor,iiu --needtodelta --needtoskiplist --zipfian1 
#./codecs --codecs pfor,iiu --needtodelta --needtoskiplist --uniformsparseclassic # ???
#./codecs --codecs pfor,iiu --needtodelta --needtoskiplist --uniformdenseclassic
#./codecs --codecs pfor,iiu --needtodelta --needtoskiplist --clustersparseclassic # ???
#./codecs --codecs pfor,iiu --needtodelta --needtoskiplist --clusterdenseclassic
#./codecs --codecs pfor,iiu --needtodelta --needtoskiplist --uniformsparse
#./codecs --codecs pfor,iiu --needtodelta --needtoskiplist --uniformdense
#./codecs --codecs pfor,iiu --needtodelta --needtoskiplist --clustersparse
#./codecs --codecs pfor,iiu --needtodelta --needtoskiplist --clusterdense
./codecs --codecs pfor,iiu --needtodelta --needtoskiplist --zipfian1
#./codecs --codecs pfor,iiu --needtodelta --needtoskiplist --zipfian2
#./codecs --codecs pfor,iiu --needtodelta --needtoskiplist --vclusterdynamic # XXX
#./codecs --codecs pfor,iiu --needtodelta --needtoskiplist --crazyclusterdynamic # XXX
#./codecs --codecs pfor,iiu --needtodelta --needtoskiplist --clusterdynamicsmall
#./codecs --codecs pfor,iiu --needtodelta --needtoskiplist --uniformdynamicsmall
#./codecs --codecs pfor,iiu --needtodelta --needtoskiplist --clusterdynamic
#./codecs --codecs pfor,iiu --needtodelta --needtoskiplist --uniformdynamic
#./codecs --codecs pfor,iiu --needtodelta --needtoskiplist --clusterdynamicpredelta
#./codecs --codecs pfor,iiu --needtodelta --needtoskiplist --uniformdynamicpredelta
#./codecs --codecs pfor,iiu --needtodelta --needtoskiplist --sillyuniformdynamic # XXX

