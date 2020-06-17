# SARSCoV2_phylodynamics

The bash and R scripts in this repo can be used for the following tasks:
1) Downsample a global SARS-CoV-2 FASTA and then further downsample sequences exogenous to a region of interest in a time-stratified manner. 
1) Randomly select exogenous sequences, including the closest cophenetic neighbors to sequences from the region of interest.
1) Summarize and visualize cumulative incidence, prevalance, and instantaneous incidence from BEAST2::Phydyn trajectory files.

The Exog_seq_selection bash script is a pipeline that will call the other scripts. A FASTA file of global sequences (global.fasta) needs to be supplied in the output directory. The output directory and a code for the region of interest (e.g. 'Dane' or 'Milwaukee') also need to be piped into the bash script. 

\*\*Full credit for filter_fasta_by_list_of_headers.py goes to StackExchange user 
[Kamil S Jaron](https://bioinformatics.stackexchange.com/users/57/kamil-s-jaron)
who provided the script in the StackExchange post found [here](https://bioinformatics.stackexchange.com/questions/3931/remove-delete-sequences-by-id-from-multifasta).
