# JCV
JCV is an R script for producing JCV matrix, heatmap from list of (species,COG) data file as well as noa and sif files for Cytoscape.

The JaccardClusters.R script reads a list of (species, COG) data pairs and calculates the Jaccard Coefficient Value (JCV) between all possible species pairs.

The JCV is defined as the intersect divided by the union of the gene content between two species:

The JCV is calculated in the following way: 

JCV = |AnB|/(|A|+|B|-|AnB|)

that is, the intersection of common genes divided by the union of all genes for species A and B, where 0 <= J <= 1.

A sample input looks like the following:

species                                 orthologue ID
Acanthamoeba_castellanii_mamavirus      NCVOG0001
Acanthamoeba_castellanii_mamavirus      NCVOG0001
Acanthamoeba_castellanii_mamavirus      NCVOG0001
Acanthamoeba_polyphaga_mimivirus        NCVOG0001
Acanthamoeba_polyphaga_mimivirus        NCVOG0001
Acanthamoeba_polyphaga_mimivirus        NCVOG0001
Acanthamoeba_polyphaga_moumouvirus      NCVOG0001
Acanthamoeba_polyphaga_moumouvirus      NCVOG0001
Bathycoccus_sp._RCC1105_virus_BpV1      NCVOG0001
Cafeteria_roenbergensis_virus_BV-PW1    NCVOG0001

The script can be run in the following way:

Rscript JaccardClusters.R input.txt

The script produces four output files:

jcv_clusters.mx - the JCV matrix for all species pairs, which is used to make the heatmap
jcv_clusters.jpg - the JCV heatmap which graphically depicts JCVs for all species pairs. Lighter colors represent JCV closer to 1.0 (continuity), darker ones represent JCVs closer to 0.0 (discontinuity)
jcv_clusters.sif - file useful for Cytoscape
jcv_clusters.sif - file useful for Cytoscape
