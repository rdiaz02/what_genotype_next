#!/bin/bash


## We are careful spacing some runs to prevent huge loads.

## hard to tell which is the best number for CORES. 20 seems good.
CORES=20
export OMP_NUM_THREADS=$CORES
export OPENBLAS_NUM_THREADS=$CORES

## Do not launch all at the same time. MHN seems to use up to 5 threads
## and it can finish in about 1200
for DATASET in {1..9}
do
    nohup R-devel --vanilla --slave -f run-all-biol-examples-max-genes-15.R --args $CORES $DATASET &> run-all-biol-examples-$CORES-$DATASET-max-genes-15.Rout &
    sleep 60
done

for DATASET in {10..14}
do
    nohup R-devel --vanilla --slave -f run-all-biol-examples-max-genes-15.R --args $CORES $DATASET &> run-all-biol-examples-$CORES-$DATASET-max-genes-15.Rout &
    sleep 240
done

for DATASET in {15..17}
do
    nohup R-devel --vanilla --slave -f run-all-biol-examples-max-genes-15.R --args $CORES $DATASET &> run-all-biol-examples-$CORES-$DATASET-max-genes-15.Rout &
    sleep 1000
done

DATASET=18
nohup R-devel --vanilla --slave -f run-all-biol-examples-max-genes-15.R --args $CORES $DATASET &> run-all-biol-examples-$CORES-$DATASET-max-genes-15.Rout

## Make sure all smaller data sets done
for DATASET in {19..27}
do
    nohup R-devel --vanilla --slave -f run-all-biol-examples-max-genes-15.R --args $CORES $DATASET &> run-all-biol-examples-$CORES-$DATASET-max-genes-15.Rout &
    sleep 1200
done
