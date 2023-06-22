#!/bin/bash
input="${1:-/dev/stdin}"
trees=$2
# ids=$3

> $trees
#> alns_121.txt

while read line; do

duplicates=$(echo $line | cut -f4 -d' ' | nw_labels -I - | cut -f2 -d'_' | sort | uniq -d)

if [ -z "$duplicates" ]; then
	echo $line | cut -f4 -d' '  >> $trees
	# echo $line | cut -f1 -d' '  >> $ids
fi

done < $input
