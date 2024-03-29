# This script takes a fasta file as input and searches for the most parsimonious
# tree using random one random stepwise addition followed by NNI branch swaping.
# It expects the phylommand programs to be found in the PATH.
# The script is a prof of concept rather than an efficient way of finding the MP
# tree.

TREES=($(treeator -r -v -d $1))
BESTSCORE=$(echo ${TREES[0]} | treeator -p -d $1 | sed 's/^\[score: //' | sed 's/\].*//')
PREV_BEST=$((BESTSCORE+1))
TREES_IN_MEM=100
NNI_RUN=1
echo "Stepwise addition score: " $BESTSCORE
while [ $BESTSCORE -lt $PREV_BEST ]
do
    PREV_BEST=$BESTSCORE
    echo "Doing NNI on ${#TREES[@]} trees"
    NNI_RUN=$((NNI_RUN+1))
    for j in ${!TREES[*]}
    do
	NNITREES=( $(echo ${TREES[$j]} | treebender -0 -v --nni all) )
	echo "NNI on tree $j (${#NNITREES[@]} trees)"
	for i in ${!NNITREES[*]}
	do
	    #echo 'Test tree'
	    SCORE=$(echo ${NNITREES[$i]} | treeator -p -d $1 | sed 's/^\[score: //' | sed 's/\].*//') 
	    #echo $i
	    #echo $SCORE
	    if [ $SCORE -lt $BESTSCORE ]
	    then
		BESTSCORE=$SCORE
		unset BESTTREES
		BESTTREES=( ${NNITREES[$i]} )
		echo "Found better tree, saving only that tree"
		echo "Updated best score: " $BESTSCORE
	    elif [ $SCORE -eq $BESTSCORE ] && [ ${#BESTTREES[@]} -lt $TREES_IN_MEM ]
	    then
		N_EQUAL=0;
		for k in ${!BESTTREES[*]}
		do
		    if [ $(echo ${NNITREES[$i]} ${BESTTREES[$k]} | contree -r | head -1 | sed "s/^.*: //" | sed "s/ (.*)//") -eq 0 ]
		    then
			N_EQUAL=$((N_EQUAL+1))
		    fi
		done
		if [ $N_EQUAL -eq 0 ]
		then
		    BESTTREES+=( ${NNITREES[$i]} )
		    echo "Equally parsimonious tree found (${#BESTTREES[@]})"
		fi
	    fi
	done
    done
    TREES=( ${BESTTREES[@]} )
    echo "Number of most parsimonious trees: ${#TREES[@]}"
done
for i in ${!TREES[*]}
do
    echo ${TREES[$i]}
done
echo $BESTSCORE
