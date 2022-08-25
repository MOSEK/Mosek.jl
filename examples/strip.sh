#!/bin/bash

cd $(dirname $0)

for f in *.jl; do
	sed '/TAG:/d' < $f \
		| sed "s/\$\${file}/$f/g" \
		| sed 's/\$\$copyright/Copyright (c) MOSEK ApS, Denmark. All rights reserved./g' \
		> tmp
	mv tmp $f
done	
