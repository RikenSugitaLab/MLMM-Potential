#!/bin/bash
for i in `seq 1 21`;
  do echo $i
  sed -e "s/ID/${i}/g" submit_read_data.sh > aa_${i}.sh
  qsub aa_${i}.sh
done

