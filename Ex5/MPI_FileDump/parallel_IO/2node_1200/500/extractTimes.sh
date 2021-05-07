#!/bin/bash


for procs in 2 4 8 16 20
 do
	grep -r 'total elapsed time on processor' ${procs}_1200_hdf5.txt > ${procs}_hdf5_totTime.txt
	grep -r 'communication elapsed time on processor' ${procs}_1200_hdf5.txt > ${procs}_hdf5_commTime.txt
        grep -r 'computation time for the evolving the bulk by processor' ${procs}_1200_hdf5.txt > ${procs}_hdf5_calcTime.txt
        grep -r 'computation time for dumping the files by processor' ${procs}_1200_hdf5.txt > ${procs}_hdf5_dumpTime.txt
 done


mkdir commTime
mkdir totTime
mkdir calcTime
mkdir dumpTime

mv *commTime.txt commTime/
mv *calcTime.txt calcTime/
mv *totTime.txt totTime/
mv *dumpTime.txt dumpTime/
