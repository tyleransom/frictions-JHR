#!/bin/sh

for item in cflFrict25*.m

do
	filename=$(basename "$item")
	item_stem="${filename%.*}"
	echo "$item_stem"
	mathup3 ${item_stem}.m ${item_stem}.log
	# mathup2 ${item_stem}.m ${item_stem}.log
	# if [[ $(($i%1)) -eq 0 ]]
	# then
	# 	sleep 75m;
	# fi
	# i=$((i+1));
done
