#!/bin/bash


for procs in 2 4 8 16 20
 do
	grep -r 'total elapsed time on processor' NB_${procs}.txt > NB_${procs}_totTime.txt
	grep -r 'communication elapsed time on processor' NB_${procs}.txt > NB_${procs}_commTime.txt
        grep -r 'computation time for the evolving the bulk by processor' NB_${procs}.txt > NB_${procs}_calcTime.txt
 done


mkdir commTime
mkdir totTime
mkdir calcTime


mv *commTime.txt commTime/
mv *calcTime.txt calcTime/
mv *totTime.txt totTime/
