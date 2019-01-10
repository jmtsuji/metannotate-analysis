# metannotate-analysis
Collection of personal scripts to process output from [MetAnnotate (DoxeyLab tool)](https://bitbucket.org/doxeylab/metannotate)

## 1. metannotate_barplots.R
Library of functions to generate bar charts or bubble charts of major taxa in a sample based on taxonomic or functional genes.
**This is an early version of this script (beta) -- don't trust this yet for real analyses without contacting me about usage. I'd appreciate any feedback on potential issues or feature requests.**

### Dependencies (R packages)
```
library(argparser)
library(parallel)
library(futile.logger)
library(roxygen2)
library(tools)
library(glue)
library(tibble)
library(plyr)
library(dplyr)
library(ggplot2)
library(reshape2)
library(scales)
```
Install R packages within R via, `install.packages("argparser")`, for example. Alternatively, you could also set up a conda environment, for those familiar with conda.

### Usage
Source the script from your R console to then be able to use the included functions to analyze your data in an explorative manner:
```
# **Link may not yet be correct!
source(https://github.com/jmtsuji/metannotate-analysis/archive/v0.9.5.metannotate_barplots.R)
```

In addition, before you get started, you need to make sure that you have the following files/information from MetAnnotate:
- all_annotations_[etc].tsv -- this can be either gzipped or not
- lengths of all HMMs used -- you can get this from the header of the HMM if you look at the first few lines in a text editor.

From here, you can then analyze your data a number of different ways. The recommended workflow is below:

#### 1. Load the data
```
metannotate_table_filename <- "path_to_your_all_annotation_tsv_file"
metannotate_data <- read_metannotate_tibble(metannotate_table_filename)
```

#### 2. Print setup templates
```
setup_templates <- create_setup_templates(metannotate_data, write_tables = TRUE)
```
This will write two tables to your working directory: `hmm_info_template.tsv` and `dataset_info_template.tsv`
You need to open these (e.g., in Excel) and fill out the required info to guide subsequent data processing and plotting.

`hmm_info_template.tsv` will look something like:
```
raw_name	HMM.Family	HMM_length	notes
TIGR01115_9			
TIGR01281_10			
TIGR02019_11			
cyc2_0			
dsrA_4			
pmoA_12			
soxB_14			
BChl_A_8			
nifH_new_7			
rpoB_13			
mcrA_5			
mmoX_6			
```

You need to fill in:
- `HMM.Family`: The readable name you want to give the HMM (e.g., rpoB)
- `HMM_length`: The length of the HMM, available in the HMM header (for normalization calculations)
- The _order_ of the HMMs in this table dictates their order in the output plot
- You can optionally _omit_ HMMs you aren't interested in from the table. They will be dropped during data processing.
- `notes` is just for your own reference if you want to write something down.

Example filled out table (note how I also moved around and deleted some rows):
```
raw_name	HMM.Family	HMM_length	notes
rpoB_13	rpoB	2842	This is the taxonomic marker gene I plan to use
cyc2_0	cyc2_3GSB	411	
dsrA_4	dsrA	369	
soxB_14	soxB	570	
BChl_A_8	fmoA	369	
TIGR01115_9	pufM	570	
TIGR01281_10	bchL	369	
TIGR02019_11	bchJ	570	
```

`dataset_info_template.tsv` will look something like:
```
raw_name	Dataset
L227_2014_6m_QC_R1_frag_2	
L227_2013_6m_QC_R1_frag_0	
L442_2011_16_5m_QC_R1_frag_6	
L227_2013_8m_QC_R1_frag_1	
L227_2014_8m_QC_R1_frag_3	
L442_2014_15m_QC_R1_frag_7	
L227_S_6D_QC_R1_frag_4	
L304_S_6D_QC_R1_frag_5	
```

You need to fill in:
- `Dataset`: The readable name you want to give the HMM (e.g., rpoB)
- The _order_ of the datasets in this table dictates their order in the output plot
- You can optionally _omit_ datasets you aren't interested in from the table. They will be dropped during data processing.

Example filled out table (note how I also moved around and deleted some rows, and even used some "special characters" in the Dataset names):

```
raw_name	Dataset
L227_2013_6m_QC_R1_frag_0	L227 2013 6m
L227_2013_8m_QC_R1_frag_1	L227 2013 8m
L227_2014_6m_QC_R1_frag_2	L227 2014 6m
L227_2014_8m_QC_R1_frag_3	L227 2014 8m
L442_2011_16_5m_QC_R1_frag_6	L442 2011 16.5m
L442_2014_15m_QC_R1_frag_7	L442 2014 15m
```

Finally, once done, back in R, map the information from these tables onto your dataset:
```
hmm_naming_info_filename <- "path_to_your_hmm_naming_file"
dataset_naming_info_filename <- "path_to_your_dataset_naming_file"
metannotate_data_mapped <- map_naming_information(metannotate_data, hmm_naming_info_filename, dataset_naming_info_filename)
```

Now, you're ready to go for downstream analysis

#### 3. Iteratively explore the MetAnnotate data
Use the convenient wrapper function `explore_metannotate_data` to try a number of e-value cutoffs, thresholds for plotting taxa, and so on, to explore your data.
This function can ultimately be used to make nearly publication-ready plots using the `colouring_template_filename` described later.

Example usage:
```
metannotate_plot <- explore_metannotate_data(metannotate_data_mapped, evalue = 1e-10, taxon = "Family",
                         normalizing_HMM = "rpoB", plot_type = "bar", 
                         top_x = 10, percent_mode = "within_sample",
                         colouring_template_filename = NA)
print(metannotate_plot)
```

What this does:
1. _Filters_ your data to an e-value threshold of your choice (defined in `evalue`). Will then print helpful summary tables to the screen to show the % change in the # of hits before and after e-value filtration. Future versions of this script may also help to provide guidance on the selection of an appropriate e-value cutoff, but for now, the user must decide on this themselves.
2. _Collapses_ the data to a given taxonomic rank (MUST be one of: domain, phylum, class, order, family, genus, species; case insensitive), as defined in `taxon`. This sums up hit counts to each taxon at the given rank, creating something kind of like an 'OTU table' in 16S amplicon analysis.
3. A: _Normalizes_ the data by HMM length. Longer HMMs tend to have more hits than shorter HMMs from raw read data (with generally linear correlation between HMM length and hit bias), so this script divides the number of HMM hits by the length of the HMM. This allows hit counts from different HMMs to be cross-compared. You already provided the HMM lengths in the `hmm_info_template.tsv` file described above.
3. B: _Normalizes_ the data according to the total sequencing depth of each sample. This allows for HMM hit counts to be meaningfully compared between metagenome datasets. Total sequencing depth is defined by the total (lenth-normalized) hits to a single-copy taxonomic marker gene (defined in `normalizing_HMM`). As such, your MetAnnotate table MUST contain a HMM of a single-copy taxonomic marker gene like __rpoB__ or __dnaK__. `normalizing_HMM` should be the "HMM.Family" name that you gave to your normalizing HMM in the `hmm_info_template.tsv` file described above.
As a consequence of this double-normalization, the output data for each functional gene is expressed as the relative abundance of the gene compared to the single copy taxonomic marker gene. For example, double normalization might show that total __dsrA__ hits represent 20% of total __rpoB__ hits, leading you to conclude that __dsrA__ must be a fairly commonly held gene in the whole microbial community in your sample. You might see that the __dsrA__ hits classified to the __Chlorobiaceae__ family represent 10% of total __rpoB__ hits (or 50% of total __dsrA__ hits), leading you to conclude that __Chlorobiaceae__ are prominant sulfur-cycling microorganisms in the system.
However, an important caveat: HMM length normalization attempts to allow HMM hit values to be directly cross-compared, but it CANNOT account for the inherent bias of different HMMs. A stringent HMM will still get fewer hits than a relaxed HMM, based on the probability frequencies defined in the HMM profile. As such, you should not 'hang your hat' on between-HMM comparisons that this script outputs. Within HMM comparisons are likely reliable, but between-HMM comparisons may be biased. So be careful before saying that one gene is more/less prevalent than another if the #s are close.
Once normalization is finished, the function prints some normalization stats to the screen for the user's interest.
4. _Subsets_ the data to the most abundant taxa, for plotting purposes.
- `top_x`: if >=1 (e.g., 10), then the script subsets the top __ most abundant taxa within each sample for plotting. If <1 (e.g., 0.01), then the script subsets all taxa of __ proportional abundance or higher within each sample for plotting.
- `percent_mode`: If `top_x` <1 (i.e., in proportional abundance mode), then there are two different ways to subset by proportional abundance for functional genes. Specify the preferred method here. If `within_sample` is selected, then the script will subset all taxa with __ proportional abundance or higher _based on the proportional abundance of that taxon in the normalizing_HMM data_. If `within_HMM` is selected, then the script will subset all taxa with __ proportional abundance or higher _based on the proportional abundance of the taxa within each HMM_. The main case where `within_HMM` is helpful is when one functional gene HMM accounts for a very small proportion of the total hits in the dataset, but you still want to see what taxa are there. Play around with these settings yourself to test them out.
5. _Plots_ the data. Can do this either as a "bar" plot or a "bubble" plot, as defined in `plot_type`. If you want to make the plot look more beautiful, you can use the `colouring_template_filename` feature described below.

Play around with the script parameters until you are satisfied. Then, when you want to make a finalized plot, improve plot colours and save to a PDF as described below.

#### 4. Beatify and export the plot
Once you have settings for `explore_metannotate_data` that you are satisfied with, you can change the plot colours to be more meaningful.

Run your plot's code again, but specify a save location for the `colouring_template_filename` instead of `NA`. This file should NOT already exist. E.g.,
```
metannotate_plot <- explore_metannotate_data(metannotate_data_mapped, evalue = 1e-10, taxon = "Family",
                         normalizing_HMM = "rpoB", plot_type = "bar", 
                         top_x = 10, percent_mode = "within_sample",
                         colouring_template_filename = "colouring_template.tsv")
```

This will write `colouring_template.tsv` to your working directory. The file looks something like (if `order` was the subsetting taxon):
```
order	colour
Methanomicrobiales	#F57A5F
Methanotrichales	#E88521
Holophagales	#D79000
Bacteroidales	#C29A00
Chlorobiales	#A9A400
Chloroflexales	#89AC00
Desulfobacterales	#5DB300
Elusimicrobiales	#00B826
Myxococcales	#00BD60
Phycisphaerales	#00C086
Azospirillales	#00C1A7
Betaproteobacteriales	#00BFC4
Methylococcales	#00BBDD
Steroidobacterales	#00B3F2
RFP12	#00A7FF
Chthoniobacterales	#7299FF
Pedosphaerales	#AB88FF
Syntrophales	#D177FF
```

Auto-generated HTML colour codes are provided, as well as an auto-generated sort order. You can now change these to HTML colour codes of your choise, and you can modify the order of the rows to plot the taxa in the order you prefer. For example:
```
order	colour
Chlorobiales	#339933
Chloroflexales	#89AC00
Desulfobacterales	#5DB300
Methylococcales	#00BBDD
Methanomicrobiales	#F57A5F
Methanotrichales	#E88521
Holophagales	#D79000
Bacteroidales	#C29A00
Elusimicrobiales	#00B826
Myxococcales	#00BD60
Phycisphaerales	#00C086
Azospirillales	#00C1A7
Betaproteobacteriales	#00BFC4
Steroidobacterales	#00B3F2
RFP12	#00A7FF
Chthoniobacterales	#7299FF
Pedosphaerales	#AB88FF
Syntrophales	#D177FF
```

Once done, you can specify the final filename as `colouring_template_filename` and then run again. For example, if you renamed your final file `colouring_guide.tsv`:
```
metannotate_plot <- explore_metannotate_data(metannotate_data_mapped, evalue = 1e-10, taxon = "Family",
                         normalizing_HMM = "rpoB", plot_type = "bar", 
                         top_x = 10, percent_mode = "within_sample",
                         colouring_template_filename = "colouring_guide.tsv")
print(metannotate_plot)
```

Once you are satisfied with the plot, you can save it using:
```
output_filename <- "path_to_output_file.pdf"
plot_width <- 200
plot_height <- 300
ggsave(file = output_filename, width = plot_width, 
         height = plot_height, units = "mm")
```
Play around with the width and height to get it to your liking.

Done! You can do further fine-scale edits in a program like Inkscape.


### How to interpret the plot (and some nitty gritty details)
Two normalization steps are performed during the production of the plot:
1. Normalize by HMM length: longer HMMs get more hits than shorter ones (e.g., due to it overlapping a larger proportion of a genome and so hitting more short reads). Thus, this script divides hit totals by HMM length (assuming a linear relationship between length and hit numbers) to attempt to account for this bias. This allows for comparison between HMMs within a single metagenome.
2. Normalize by total marker gene (e.g., _rpoB_) hits within each sample: each metagenome will have a slightly different number of relevant reads. To account for this difference, one can express all HMM hits within a single metagenome as relative abundances to a single-copy taxonomic marker gene that has predictable behaviour between different environments. This script sums the total number of hits for the given taxonomic marker gene within each sample (AFTER length normalized) and then divides other length-normalized HMM hits by this number in order to express them as proportional abundances relative to the marker gene.

This lays the framework for understanding the bar charts. For each plotted HMM:
* Total hits relative to the taxonomic marker are shown as a grey bar ("#808080"). This gives an indication of the **total abundance of that gene within the microbial community relative to the marker**. For example, if the grey bar is at 30% for the _nifH_ gene, then potentially, ~30% of microorganisms within the community possess that gene (or 15% possess two copies of that gene, and so on).
* The top specified taxa are shown as coloured bars (as specified by the user)

### Future development plans
- Make runnable via command line
- Allow for a more carefully fine-scale user workflow for advanced users. For example, make a way to export intermediate tables or the final data table, not just the plot.
- Turn this into a proper R package?

