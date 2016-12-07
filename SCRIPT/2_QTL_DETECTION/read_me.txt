Bonjour,

Ici je mets le scripts de détection QTL.
Il faut le commenter en anglais, le mettre plus propre, et donner une ligne de code example.

Calcul pour DS:
qsub -q normal.q -b yes -cwd -N tmp_qtl "Rscript /NAS/g2pop/HOLTZ_YAN_DATA/programmes/QTL_detection_with_QTLRel.R  genotypage.csv phenotypage.csv  carte"



