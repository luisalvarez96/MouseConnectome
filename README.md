# Codes to reproduce analysis of Mouse Connectome project. in preparation

Overall idea of the project consists in utilizing the observed cluster synchronization in the Mouse brain, thorugh cfos data, to understand the imnpact of the connectome structure in the learning/memory process of Mouse. 
We seek to identify a Region(s) of Interest (ROIs) whose removal from the connectome will have a measuring impact on the Mouse's learning/memory capacity, without these ROIs being obvious.
The structural connectome is reconstructed from the Allen connectome data. 
We then use a symmetry repair algorithm, utilizing gurobi to perform optimization, to reconstruct the functional connectome based on the observed cluster synchronization from cfos experimental data.

**cleanish.R** 
