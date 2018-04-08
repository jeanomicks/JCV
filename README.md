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

=====Version 3 of Jaccard Coefficient script (April 7, 2018)=====

The latest version of the R script does several new things:
- it creates an output directory for the output files (for this you have to add an output directory name as the fourth parameter)
- it creates the heatmap in red to yellow color
- it calculates those genes which belong to the core genome and the pan genome of the cluster
   o the core genome being the collection of genes common to all species in the cluster
   o the pan genome being the collection of genes in at least one species of the cluster
- it adds the sizes of the core and pan genome and the core/pan genome ratio to the stats file

=====Multi-algorithm Jaccard Coefficient Method=====

The original version of the Jaccard Coefficient Method predicts clusters/baramins based on k-means clustering.
Two extra clustering algorithms are available to choose from when using the multi-algorithm Jaccard Coefficient Method.
For this use the JaccardCoefficientMulti.R script. These two extra algorithms determine clusters with at least three members.

These 2 algorithms are:

1.) the PGQ (pgq) method: here the algorithm determines clusters based on a high pan-genome quotient (PGQ). A PGQ can be detected when the algorithm starts out from a seed species and keeps adding newer and newer species. The PGQ is tracked constantly until there is a sharp drop in the value, denoting that all members have been added to the cluster which highly overlap in gene content with each other.

test run:
Rscript JaccardClustersMulti3.R pgq input2.txt outlier

2.) The Matrix Cut (mxcut) method: Here the user has to supply a JCV cutoff value. The algorithm determines members of a cluster whose members each have a JCV with each other in a pairwise manner. In other words, it takes the whole JCV matrix in graph representation and removes edges (interspecies relationships) which have a JCV below the cutoff value. The remaining "cut" graph consists of the predicted clusters.

test run: 
Rscript JaccardClustersMulti3.R mxcut input2.txt outlier 0.7
