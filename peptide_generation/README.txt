To run the peptide generation for a particular sample, execute

python2 -u ./compute_peptide_lists.py --parallel_rank <RANK> --run_id <ID_OF_THE_RUN>

The parallel rank RANK is a running integer in the interval [0,number of samples to process) where
<number of samples to process> is the length of the list that was specified as a DONOR_PATH
to the script, which enables easy parallel processing in a cluster environment. The run id
ID_OF_THE_RUN defines the folder name of the output peptides. Use the command line
flag '--batch_mode' to re-direct stdout/stderr to log files in a cluster environment.
Use '--gtex_normals' to generate GTEX background outputs. GRAPH_DATA and GRAPH_METADATA_PATH
then also need to be modified accordingly to pass the correct splice graphs for GTEX
peptide generation mode.

The script <compute_peptide_lists.py> needs to be adapted to your local environment by
substituting paths for the input/output data and the location of the spladder dependency
at the bottom of the script. Default paths will be generated from this. Alternatively,
the location of input/output paths can be passed directly as command line arguments to the script.

For generating the foreground/background peptides, cluster wrapper scripts
{generate_peptides_gtex_normals.sh, generate_peptides_tcga.sh} were used and are
provided with the code.





