## Data source
Ebolavirus and marburgvirus multialignment file source: 
http://hgdownload.soe.ucsc.edu/goldenPath/eboVir3/multiz160way/

## File preparation
This multialignment contains a subset - only one block - of the original multialignment. This block is long enough to give interesting results when used as PangTreeBuild input and is short enough to be visualised by PangTreeVis smoothly.
In order to be easily processable by PangtreeBuild, also the sequences coordinates were changed. 

Also, sequences coordnates were changed in order to avoid the need of missing nucleotides refillment. It is imitated that this multialignment contains full sequences.

## Recommended PangtreeBuild parameters for this file
- multialignment: "example_data/Ebola_subset/ebola_subset_774_indexes.maf",
- metadata: "example_data/Ebola_subset/metadata.csv",
- affinity: "tree",
- output_full