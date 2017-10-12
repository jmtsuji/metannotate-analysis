# metannotate-analysis
Collection of personal scripts to process output from [MetAnnotate (DoxeyLab tool)](https://bitbucket.org/doxeylab/metannotate)

## 1. metannotate_barplots.R
Script to generate bar charts of major taxa in a given sample with a given functional gene.
**This is an early version of this script (alpha) -- don't trust this yet for real analyses without contacting me about usage. I'd appreciate any feedback on potential issues or feature requests.**
### Usage
This script must be run in several steps from within a console (e.g., RStudio), modifying paramters in the "User variables" section in the first few lines of the code. The input is the "annotations" TSV file provided as output from MetAnnotate. This file can contain data for multiple HMMs and samples; the script will plot all (or a subset, if you desire) in the same figure.
1. After selecting your working directory and MetAnnotate summary table, run ```script_setup <- TRUE```. This will print out two tables to fill in information about your samples and HMMs needed for further processing.
* ```hmm_info_template.tsv``` - Provide the actual name of the gene for each HMM (will show up in the final plot; white spaces okay). Also, provide the amino acid length of each HMM. This information is available in the header of each HMM (e.g., if opened in a ext editor; the LENG parameter) and is critical for normalizing the data.
* ```sample_info_template.tsv``` - Provide the actual name of each sample (i.e., "Dataset" in MetAnnotate terminology) that will appear in the final figure (white spaces are okay). **The order of the samples in this file will dictate their order in the final plot.** Also, if you ran unmerged paired-end reads through MetAnnotate, and now have "raw" sample names for both the forward and reverse reads, you can combine them (if you'd like) by giving both the forward and reverse reads the same "final" name.

Once finished, set the filled-in template files as inputs in the "Required supplemental data" section.

2. Now, you can explore your data by tweaking the "Basic plot" settings and making plots iteratively:
* ```HMMs_to_plot``` - provide a subset of HMMs that you'd like to appear on the final figure (use the final plotting names listed in the ```hmm_info_template.tsv``` file). **The order in which you provide the HMM names dictates their order in the final plot.**
* ```normalizing_HMM``` - this is important. In order to compare MetAnnotate gene hit results between metagenomes ("Datasets"), you must normalize gene hits within each metagenome to some sort of **single-copy taxonomic marker gene (e.g., _rpoB_)**. In this way, you represent the hits of each gene as relative abundances within a sample compared to the abundance of the marker gene. Specify the name of that marker gene here (HMM name needs to match that provided in the ```hmm_info_template.tsv``` file).
* ```tax_rank_to_plot``` - e.g., Superkingdom, Phylum, Class, Order, Family, Genus, or Species. Script will sum together all hits at that given taxonomic rank and then plot as stacked bars.
* ```top_number_to_plot``` - It is too cluttered to include all MetAnnotate hits in a bar chart. This script subsets those hits to keep only the most abundant taxa. Here, you decide what the criteria are for choosing those taxa. If < 1, the script will plot all taxa with at least this relative abundance in the community. (E.g., if 0.01, will plot all taxa with >= 1% relative abundance to the taxonomic marker gene.) If > 1 (e.g., 10), then plot the top ___ (e.g., 10) most abundant taxa for each gene. One tends to need to optimize this for each dataset to get an informative yet not overwhelming number of taxa to display.

**This script will generate barplots with automatic colours and taxa orders. Once you find a layout you really like, you can then touch up the plot using the advanced settings (4).**

3. Specify some basic output parameters in the "Options for output" section, such as whether or not to print summary tables from key steps, or whether or not to print PDFs.

4. Making (nearly) publication-ready custom plots
Once you have found a set of "basic" plotting settings that you are satisfied with, you can overlay more custom user information to make the plot more meaningful and visually appealing. This must be done in a series of steps:
* a. Set ```print_custom_plot_template <- TRUE```. This will cause the script to print off a template file with the current "basic" plot settings for the user to fill out.
* b. Fill out the template file. Sort the taxa (i.e., rows) in the order that you want them to appear in the plot. It helps to group taxa with similar (potential) functional roles together, for example. Also, specify the colour that you want each taxon to be plotted with. You could potentially group taxa by potential function (or phylogenetic similarity) and put into groups of related colours, for example. As a related example, see [this paper, Figure 3](http://doi.org/10.1038/srep46708).
* c. Once done, specify the completed template file as ```custom_plot_template_filename``` and set ```print_custom_plot_template <- FALSE``` again. Now, set ```print_custom_plot <- TRUE``` to generate your custom plot.
* d. You'll probably still need to tweak this plot a bit after it is produced (e.g., using a PDF editor like [Inkscape](https://inkscape.org/en/). Note that you can also dig into the script and change the plot dimensions and so on, if you'd like.

### How to interpret the plot (and some nitty gritty details)
Two normalization steps are performed during the production of the plot:
1. Normalize by HMM length: longer HMMs get more hits than shorter ones (e.g., due to it overlapping a larger proportion of a genome and so hitting more short reads). Thus, this script divides hit totals by HMM length (assuming a linear relationship between length and hit numbers) to attempt to account for this bias. This allows for comparison between HMMs within a single metagenome.
2. Normalize by total marker gene (e.g., _rpoB_) hits within each sample: each metagenome will have a slightly different number of relevant reads. To account for this difference, one can express all HMM hits within a single metagenome as relative abundances to a single-copy taxonomic marker gene that has predictable behaviour between different environments. This script sums the total number of hits for the given taxonomic marker gene within each sample (AFTER length normalized) and then divides other length-normalized HMM hits by this number in order to express them as proportional abundances relative to the marker gene.

This lays the framework for understanding the bar charts. For each plotted HMM:
* Total hits relative to the taxonomic marker are shown as a grey bar ("#808080"). This gives an indication of the **total abundance of that gene within the microbial community relative to the marker**. For example, if the grey bar is at 30% for the _nifH_ gene, then potentially, ~30% of microorganisms within the community possess that gene (or 15% possess two copies of that gene, and so on).
* The top specified taxa are shown as coloured bars (as specified by the user)



