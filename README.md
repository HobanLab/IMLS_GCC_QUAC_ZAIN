<b><p><h1 style="color:red;font-size:20px;"> Project Description</b></p></h1>

This repo contains code for determining how well wild genetic diversity is captured <i>ex situ</i> for both <i>Quercus acerifolia</i> (abbreviated to QUAC in this repo), a rare and threatened North American oak species native to Arkansas and <i>Zamia integrifolia</i> (abbreviated to ZAIN in this repo), the only native cycad to the continental US, is also rare and threatened species. These species are similarly threatened because they have limited population numbers, a history of fragmented habitat due to intensive human land use, and are both very well-preserved in botanic gardens. QUAC and ZAIN both have >300 individuals <i>ex situ</i>, while only ~150 individuals <i>ex situ</i> are predicted to capture large portions of diversity in these species (Rosenberger et al., 2021), so we hypothesize: 
<ul><li>That both species' wild genetic diversity should be well-represented in botanic gardens. We can test this hypothesis by determining if:</ul></li>
<ul><ul><li>Genetic diversity in wild and garden populations does not differ significantly</li></ul></ul>
<ul><ul><li>At least 95% of all wild alleles are captured <i>ex situ</i></li></ul></ul>
<ul><ul><li>Duplicates of all wild alleles are captured <i>ex situ</i></li></ul></ul>

<ul><li>Geographic diversity of both QUAC and ZAIN is well-represented in botanic gardens, which we will determine by:</ul></li>
<ul><ul><li>Running STRUCTURE and STRUCTURE harvester on wild and garden populations to determine if all wild genetic structure is represented in botanic garden collections</li></ul></ul>
<ul><ul><li>Performing PCA on both species to determine if gardens encompass all wild genetic structure</li></ul></ul>
<ul><ul><li>Assigning botanic garden individuals to wild source populations using Geneclass 2</li></ul></ul>

We also attempt to provide recommendations for improving ex situ collections by assessing if:
<ul><li>Relatedness is higher in ex situ populations than in wild populations</ul></li>
<ul><ul><li>If relatedness is higher ex situ, this limits the ability to use these collections for restoration as offspring will be highly inbred, seriously depleting the genetic diversity and therefore usefulness of the produced offspring.</ul></ul></li>
<ul><li>Resampling in situ populations provides a better picture of diversity capture and we should edit sampling protocols in the future to better protect these species</ul></li>

The code in this repo details the analyses performed on nuclear microsatellite data for QUAC and ZAIN to test our above hypotheses and determine how we can improve guidelines for creating ex situ collections. 

<b><p><h1 style="color:red;font-size:20px;">Species Descriptions</b></p></h1>

<i>Quercus acerifolia</i> is a rare oak native to Arkansas with only four (or maybe five) wild populations and around 600-1000 wild individuals. Wild leaf tissue for this project was collected from 174 wild individuals during sampling trips by colleagues in 2019. Samples of the botanic garden individuals were collected in collaboration with 15 different botanic gardens in 2019 and resulted in 316 tissue samples for genetic data collection. Genotyping was performed with 15 nuclear microsatellite loci, with 10 being expressed sequence tag associated and 5 being genomic microsatellites. Genotyping was performed with 11 microsatellite loci for both wild and garden individuals.

<i>Zamia integrifolia</i> is the only cycad species native to the continental US and has been extirpated from parts of its range in the last century due to habitat fragmentation and overharvesting by humans. For this project, 385 individuals from 10 different botanic gardens and 751 ndividuals were collected from 25 different wild populations. Genotyping was performed with 11 microsatellite loci for both wild and garden individuals.
 
<b><p><h1 style="color:red;font-size:20px;">Folder Descriptions</b></p></h1>

This repo is divided into 2 main files - Analyses and Data_Files. 

<b>Data_Files:</b> This file contains all of the input data files for each analysis performed in the repo. There are 2 main types of data files - data frames and adegenet files. 
<ul><li><b>CSV files</b> contain the genotypes of each individual for each species in a format where individual IDs are in the first column, population names are in the second column (either wild population name or botanic garden name) and an ID for what population type the individual is (either a wild or botanic garden individual) is the third column. The fourth column onward is where the allele data are stored, with one column per allele and each consecutive column pair being the 2 alleles for one locus. Alleles are indicated with their microsatellite length.</li></ul>
<ul><li><b>Adegenet files</b> are either genepop (.gen) or arlequin files (.arp). The arlequin files are only created to be converted to the genepop files, which is done using the "arp2gen" function from diveRsity. The genepop files are imported in the R package 'adegenet' to do the majority of the genetic analyses performed in this study. The genepop files are similar to the format of the CSV files, but use commas after the individual name to denote belonging to the same population as other individuals, and then a separate line with "POP" to indicate a population break. All individuals are named with their population marker. At the top of the document all the loci used to genotype the individuals are listed- each one on a line. There is also a comment line with what populations and individuals each adegenet file. </li></ul>
<ul><li><b>Structure_Files</b>are the files used to run the program STRUCTURE (Pritchard et al. 2000). All structure files are text files, but the files deemed "_str" are the files that are initially generated from the conversion to structure function in Genalex. The files titled "_str_READY.txt" are the files STRUCTURE was actually run on. The "_str_READY.txt" file is a text file with a similar layout to the CSV files described above, just with individual names in each row and a population name assigned as numerals. There are no column headers in this text file.</li></ul>
<ul><li><b>Geneclass files</b> are the files used to run the Geneclass 2 software (Piry et al. 2004) for assignment into source populations. The input files for Geneclass are defined using "input" and are in the Genepop format as described above, but with numeric coding for alleles in a two digit format. Additionally, in this folder, the Genalex files to generate the each input Genepop file is included. Genalex was used to convert these files into the input Genepop files for the Geneclass software. The "wild" genepop files are used as the reference populations in the Geneclass software and the "garden" genepop files are used as the samples to be assigned in the software.</li></ul>

<b>Naming Conventions:</b> There are consistent naming conventions for file types in this repo, which are described brielf here. 
<ul><li>Any file named as "clean" indicates that these files were cleaned in "2_clonecheck_md_relatedness.R" script, which cleans data frames and adegenet files for individuals with too much missing data (over 25%) and clones. These data files are used for all other genetic analyses in the study. </li></ul>
<ul><li>ZAIN_og files are the original data files used in analysis, with the originally scored genetic data. ZAIN_rebinned data files have been rebinned to make results more consistent by researcher (see explanation in "rebinning analysis" script).</li></ul>
<ul><li>QUAC_wK are data files that include the wild population "Kessler Mountain" which may not be true Q. acerifolia individuals. QUAC_woK files are data files with Kessler mountain individuals removed. </li></ul>
<ul><li>allpop data files include all wild and garden populations separately, while files in the "Garden_Wild" folders with "garden_wild" in the title are where populations are separated into just their population type - garden or wild - and data files with "garden_allwildpop" is in the title indicates that all wild populations are named but all garden individuals are lumped into one population. </li></ul>
<ul><li>A file marked with "df" in the title is in the format as described above, whereas "genalex" data files are in the genalex file format. </li></ul>

Organization of Data_Files folder: 
<ul><li>Adegenet_Files</li></ul>
<ul><ul><li>Garden_Wild</li></ul></ul>
<ul><li>Data_Frames</li></ul>
<ul><ul><li>Garden_Wild</li></ul></ul>
<ul><li>Geneclass_Files</li></ul>
<ul><li>Structure_Files</li></ul>

<b>Analyses:</b> The three main folders in this folder are Analysis_RScripts, Functions, and Results. The RScripts are used to run analyses on the files in the Data_Files folder, the Functions folder contains functions created to run certain analyses, and the Results folder contains the results of those analyses. 

<b>Organization of Analyses folder:</b> 
<ul><li>Analysis_RScripts</li></ul>
<ul><ul><li>1_ZAIN_Scoring_Comparison_Barplot.R</li></ul></ul>
<ul><ul><li>2_clonecheck_md_relatedness.R</li></ul></ul>
<ul><ul><li>3_gendiv_summary_stats.R</li></ul></ul>
<ul><ul><li>4_garden_wild_comparison.R</li></ul></ul>
<ul><ul><li>5_allelic_resampling.R</li></ul></ul>
<ul><ul><li>6_IBD.R</li></ul></ul>
<ul><ul><li>7_PCA.R</li></ul></ul>
<ul><ul><li>8_Structure.R</li></ul></ul>
<ul><ul><li>9_maternal_accessions.R</li></ul></ul>
<ul><ul><li>10_geneclass_sum_stats.R</li></ul></ul>
<ul><ul><li>QUAC_clonal_propagation.R</li></ul></ul>
<ul><li>Functions</li></ul>
<ul><ul><li>accession_count.R</li></ul></ul>
<ul><ul><li>dms_degree_conversion.R</li></ul></ul>
<ul><ul><li>Fa_sample_funcs.R</li></ul></ul>
<ul><ul><li>maternal_accession.R</li></ul></ul>
<ul><ul><li>relatedness_analyses.R</li></ul></ul>
<ul><ul><li>relatedness_tests.R</li></ul></ul>
<ul><ul><li>resampling.R</li></ul></ul>

<ul><li>Results</li></ul>
<ul><ul><li>Clustering</li></ul></ul>
<ul><ul><ul><li>PCA</li></ul></ul></ul>
<ul><ul><ul><li>Structure</li></ul></ul></ul>
<ul><ul><ul><ul><li>QUAC</li></ul></ul></ul></ul>
<ul><ul><ul><ul><li>ZAIN</li></ul></ul></ul></ul>
<ul><ul><li>Garden_Wild_Comparisons</li></ul></ul>
<ul><ul><li>Scoring_Comparison</li></ul></ul>
<ul><ul><li>Sum_Stats</li></ul></ul>
