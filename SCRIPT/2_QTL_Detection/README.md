QTL Detection with QTL-REL
==========

This is the script used to do a simple marker analysis using the QTL-Rel library.
For a matter of computation time, we used this script on a calculation cluster.
To use the script on the DS population, please use this line of code:

On our server
```
qsub -q normal.q -b yes -cwd -N tmp_qtl "Rscript /NAS/g2pop/HOLTZ_YAN_DATA/programmes/QTL_detection_with_QTLRel.R  genotypage.csv phenotypage.csv  carte"
```

On a classic computer:
```
Rscript /NAS/g2pop/HOLTZ_YAN_DATA/programmes/QTL_detection_with_QTLRel.R  genotypage.csv phenotypage.csv  carte"
```



