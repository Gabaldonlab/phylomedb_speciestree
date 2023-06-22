#!/bin/bash
input="${1:-/dev/stdin}"
output="${2:-/dev/stdout}"

while read -r line; do
    seed=$(echo $line | cut -d' ' -f1)
    tree=$(echo $line | cut -d' ' -f4 | nw_labels -I - | sort | uniq | paste -s -d ',')
    echo -e "$seed\t$tree"
< $input
done > $output