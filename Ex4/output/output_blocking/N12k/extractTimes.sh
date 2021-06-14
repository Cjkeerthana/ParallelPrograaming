#!/bin/bash


for procs in 20 40 60 80 160 320
 do
	grep -r 'total elapsed time on processor' Blocking_${procs}.txt > Blocking_${procs}_totTime.txt
	grep -r 'communication elapsed time on processor' Blocking_${procs}.txt > Blocking_${procs}_commTime.txt
        grep -r 'computation time for the evolving the bulk by processor' Blocking_${procs}.txt > Blocking_${procs}_calcTime.txt
 done


mkdir commTime
mkdir totTime
mkdir calcTime


mv *commTime.txt commTime/
mv *calcTime.txt calcTime/
mv *totTime.txt totTime/
