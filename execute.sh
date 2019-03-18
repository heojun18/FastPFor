# !/bin/bash
# microbenchmark

#BENCH=(zipfian1 uniformsparseclassic uniformdenseclassic clustersparseclassic uniformsparse uniformdense clustersparse clusterdense zipfian1 zipfian2 clusterdynamicsmall uniformdynamicsmall clusterdynamic uniformdynamic clusterdynamicpredelta uniformdynamicpredelta) # not working: vclusterdynamic crazyclusterdynamic sillyuniformdynamic
BENCH=(zipfian1 zipfian2)
COMPRESSION=(pfor iiu)

make clean
make -j8
make codecs

for (( i=0; i<${#BENCH[@]}; i++ ))
do
	for (( j=0; j<${#COMPRESSION[@]}; j++ ))
	do
		echo "Run ${BENCH[$i]} (${COMPRESSION[$j]})"
		./codecs --codecs ${COMPRESSION[$j]} --needtodelta --needtoskiplist --${BENCH[$i]} 
	done
done

