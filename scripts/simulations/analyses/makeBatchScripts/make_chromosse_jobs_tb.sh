#!/bin/bash

# define template file
template1_fn="scripts/csse/templates/csse_csse_modelspec_template.rev"
template2_fn="scripts/chisse/templates/csse_chisse_modelspec_template.rev"

# make job files
for i in {101..200}; do
  echo $i
  sed "s/PLACEHOLDER/$i/g" $template1_fn > "scripts/csse/csse_csse_$i.rev"
  sed "s/PLACEHOLDER/$i/g" $template2_fn > "scripts/chisse/csse_chisse_$i.rev"
done