#!/bin/bash
for f in *.jl; do cat $f | sed -e "s/\$\$copyright/Copyright (c) 2022 MOSEK ApS/g" > tmp; mv tmp $f; done
for f in *.jl; do cat $f | sed -e "s/\$\${file}/$f/g" > tmp; mv tmp $f; done
