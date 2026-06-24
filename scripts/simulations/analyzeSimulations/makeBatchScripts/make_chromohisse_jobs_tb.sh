#!/bin/bash

# define template file
template1_fn="scripts/csse/templates/chisse_csse_modelspec_template.rev"
template2_fn="scripts/chisse/templates/chisse_chisse_modelspec_template.rev"

# make job files
for i in {1..100}; do
  echo $i
  sed "s/PLACEHOLDER/$i/g" $template1_fn > "scripts/csse/chisse_csse_$i.rev"
  sed "s/PLACEHOLDER/$i/g" $template2_fn > "scripts/chisse/chisse_chisse_$i.rev"
done