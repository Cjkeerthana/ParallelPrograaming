#!/bin/bash


for procs in 2 4 8 16 32
 do
	grep -r 'total elapsed time on processor' ${procs}_hybrid.txt > ${procs}_hybrid_totTime.txt
	grep -r 'communication elapsed time on processor' ${procs}_hybrid.txt > ${procs}_hybrid_commTime.txt
        grep -r 'computation time for the evolving the bulk by processor' ${procs}_hybrid.txt > ${procs}_hybrid_calcTime.txt
        grep -r 'computation time for dumping the files by processor' ${procs}_hybrid.txt > ${procs}_hybrid_dumpTime.txt
 done


mkdir commTime
mkdir totTime
mkdir calcTime
mkdir dumpTime

mv *commTime.txt commTime/
mv *calcTime.txt calcTime/
mv *totTime.txt totTime/
mv *dumpTime.txt dumpTime/
