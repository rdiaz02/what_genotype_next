#!/bin/bash
nohup R-devel --vanilla --slave -f fit-lmertree-specific-gamma-epist-7-w-minsize-01-Obs-infinity.R &> fit-lmertree-specific-gamma-epist-7-w-minsize-01-Obs-infinity.Rout &
nohup R-devel --vanilla --slave -f fit-lmertree-specific-gamma-epist-10-w-minsize-01-Obs-infinity.R &> fit-lmertree-specific-gamma-epist-10-w-minsize-01-Obs-infinity.Rout &
nohup R-devel --vanilla --slave -f fit-lmertree-fLandscape-rsign-peaks-recoded-7genes-w-minsize-01-Obs-infinity.R &> fit-lmertree-fLandscape-rsign-peaks-recoded-7genes-w-minsize-01-Obs-infinity.Rout &
nohup R-devel --vanilla --slave -f fit-lmertree-fLandscape-rsign-peaks-recoded-10genes-w-minsize-01-Obs-infinity.R &> fit-lmertree-fLandscape-rsign-peaks-recoded-10genes-w-minsize-01-Obs-infinity.Rout &


nohup R-devel --vanilla --slave -f fit-lmertree-specific-gamma-epist-7-w-minsize-01-Obs-infinity-sqrt.R &> fit-lmertree-specific-gamma-epist-7-w-minsize-01-Obs-infinity-sqrt.Rout &
nohup R-devel --vanilla --slave -f fit-lmertree-specific-gamma-epist-10-w-minsize-01-Obs-infinity-sqrt.R &> fit-lmertree-specific-gamma-epist-10-w-minsize-01-Obs-infinity-sqrt.Rout &
nohup R-devel --vanilla --slave -f fit-lmertree-fLandscape-rsign-peaks-recoded-7genes-w-minsize-01-Obs-infinity-sqrt.R &> fit-lmertree-fLandscape-rsign-peaks-recoded-7genes-w-minsize-01-Obs-infinity-sqrt.Rout &
nohup R-devel --vanilla --slave -f fit-lmertree-fLandscape-rsign-peaks-recoded-10genes-w-minsize-01-Obs-infinity-sqrt.R &> fit-lmertree-fLandscape-rsign-peaks-recoded-10genes-w-minsize-01-Obs-infinity-sqrt.Rout &

