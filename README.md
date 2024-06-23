# Codes to reproduce analysis of Mouse Connectome project. in preparation

Overall idea of the project consists in utilizing the observed cluster synchronization in the Mouse brain, thorugh cfos data, 
to understand the imnpact of the connectome structure in the learning/memory process of the Mouse. 
We seek to identify a Region(s) of Interest (ROIs) whose removal from the connectome will have a measuring impact on the Mouse's learning/memory 
capacity, without these ROIs being obvious.
The structural connectome is reconstructed from the Allen connectome data. 
We then use a symmetry repair algorithm, utilizing gurobi to perform optimization, to reconstruct the functional connectome based on the 
observed cluster synchronization from cfos experimental data.

**cleanish.R** contains several functions to analyze the Allen connectome data, e.g. performing percolations and collapsing/removing repeated 
ROIs divided in their sub-areas.

**fromColors.R** contains several functions to analyze and reduce the Allen connectome data, along with several functions to analyze the 
correlation data from cfos and to extract and analyze the synchronized clusters (colors). Further functions are included to analyze the resulting 
network output from the symmetry repair algorithm

**repair_direct_ar_prohibit.py** corresponds to the symmetry repair algorithm, utilizing gurobi to perform an optimization algorithm 
to satisfy that all nodes are balanced colored according to the given colors. _workpath_ and _outpath_ need to be specified, as well
as the directory list, _dirlist_, within the _workpath_ where the inputs of the codes are provided. The code takes two input files:
1) one for the network to be repaired, provided as a list of edges: source\_node, target\_node and 2) a list of the nodes with 
indicating their respective colors, i.e. the synchronized cluster to which they belong. 
