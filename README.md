# Code for *Comprehensive Analysis of Alternative Splicing Across Tumors from 8,705 Patients*
Public repository containing research code for the TCGA PanCanAtlas Splicing project accompanying
the manuscript *Comprehensive Analysis of Alternative Splicing Across Tumors from 8,705 Patients*

The code is organized in several sub-directories, which contain the following content:

**`alt_splice`**
Code concerning the detection and analysis of alternative splicing events.

    **`alt_splice/event_stats`**  
    Code for visualising the event type distributions and number of alternative splicing events detected
    per cancer type.
    
    **`complexity`**
    Detection of neo-junctions predominantly present in cancer samples but absent from annotation or the GTEx outgroup.

    **`alt_splice/outliers`**  
    Code for detection and visualization of splicing outliers.

**`neoepitopes`**  
Contains code for the visualization of results on the prediction of neo-epitopes.

**`sample_processing_rna`**  
Code for pre-processing and alignment of the RNA-Seq data, including expression counting, the
collection of usage statistics, and construction of splicing graphs.

**`sqtl`**  
Visualization scripts for the sQTL analyses.

**`tSNE_and_heatmap_visualizations`**  
Code for tSNE embeddings on alternative splicing and expression data as well as scripts 
for visualization and different highlightings.

