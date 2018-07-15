#!/bin/bash

if [[ "$#" == "0" ]]; then
	dot_files=($(ls *.dot))
else
	dot_files=("$@")
fi

png_files=()
clean_png_files() {
	echo "Cleanup of ${png_files[@]}"
	rm -f "${png_files[@]}"
}
trap clean_png_files EXIT

for i in "${!dot_files[@]}"; do
	dot="${dot_files[$i]}"
	png="${dot/%.dot/.png}"
	png_files[$i]="$png"
	dot -Tpng "$dot" > "$png" &
done
wait

feh "${png_files[@]}"
