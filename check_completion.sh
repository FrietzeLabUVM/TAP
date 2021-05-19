#!/bin/bash
#reports the root prefix for any samples that have started but not completed
for f in *start; do 
  comp=${f/start/complete}; 
  if [ ! -f $comp ]; then 
    root=${f/.start/""}; 
    echo $root; 
  fi; 
done
