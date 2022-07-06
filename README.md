<b><p><h1 style="color:red;font-size:20px;"> Project Description</b></p></h1>

This repo contains code for determining how well wild genetic diversity is captured <i>ex situ</i> for both <i>Quercus acerifolia</i> (abbreviated to QUAC in this repo), a rare and threatened North American oak species native to Arkansas and <i>Zamia integrifolia</i> (abbreviated to ZAIN in this repo), the only native cycad to the continental US, is also rare and endangered species. These species are similarly threatened because they have limited population numbers, a history of fragmented habitat due to intensive human land use, and are both very well-preserved in botanic gardens. QUAC and ZAIN both have >300 individuals <i>ex situ</i>, while only ~150 individuals <i>ex situ</i> are predicted to capture large portions of diversity in these species, which asserts both species' wild diversity should be well captured in botanic gardens. 

<i>Quercus acerifolia</i> is a rare oak native to Arkansas with only four (or maybe five) wild populations and around 600-1000 wild individuals. Wild leaf tissue for this project was collected from 174 wild individuals during sampling trips by colleagues in 2019. Samples of the botanic garden individuals were collected in collaboration with 15 different botanic gardens (Bartlett Tree Research Laboratories, Missouri Botanical Garden, The Morton Arboretum, Huntington Botanical Gardens, United States National Arboretum, Arboretum Pouyouleix, Denver Botanic Gardens, Arnold Arboretum, Trompenburg Tuinen & Arboretum, Meise Botanic Garden, Peckerwood, Zoo and BG Plze, Morris Arboretum, University of Washington Botanic Gardens, Forstbotanischer Garten Tharandt, Moore Farms Botanical Garden, Chicago Botanic Gardens) over the course of a year (2018-2019) and resulted in 316 tissue samples for genetic data collection. DNA was extracted in summer 2019 by Bailie Munoz and PCRs were run 2019 - 2020. Genetic data analyses were performed 2020 - 2021.  

<i>Zamia integrifolia</i> is the only cycad species native to the continental US and has been extirpated from parts of its range in the last century due to habitat fragmentation and overharvesting by humans. For this project, Patrick Griffith collected 387 individuals from 10 different botanic gardens (Montgomery Botanical Center, Marie Selby Botanical Gardens, The Arnold Arboretum of Harvard, The Huntington Library, Art Museum, and Botanical Gardens, Fairchild Tropical Botanic Garden, University of California at Berkeley Botanic Garden, Ganna Walska Lotusland, Harry P. Leu Gardens, Naples Botanical Garden, and Key West Tropical Forest and Botanical Garden) which were sent to the Morton Arboretum in Lisle, IL for DNA extraction and PCRs. DNA was extracted by Loren Ladd in Summer 2021 and some PCRs were performed. Emily Schumacher finished PCRs in Fall 2021, and began genetic analyses in Winter 2021. These data were compared with microsatellite data collected by Patrick Griffith and colleagues for a separate population-wide analysis of the wild genetic clustering of <i>Z. integrifolia</i> individuals (see Griffith et al., 2022). Therefore, analyses were performed to make sure scoring of microsatellites was consistent between researchers and time periods. 

The code in this repo specifically details the analyses performed on the genetic data generated from these projects (all nuclear microsatellites) and determining how well diversity is captured of these species by their large <i>ex situ</i> collections.

<b><p><h1 style="color:red;font-size:20px;">Folder Descriptions</b></p></h1>

This repo is divided into 2 main files - Analyses and Data_Files. 

<b>Data_Files:</b> This file contains all of the input data files for each analysis performed in the repo. There are 2 main types of data files - data frames and adegenet files. Data frames are CSV files which contain the genotypes of each individual for each species in a format where individual IDs are in the first column, population is in the second column (either wild population name or botanic garden name) and with a third column that has an ID for what population type the individual is (either a wild or botanic garden individual). The fourth column and on in the data frame is where the allele data are stored, with one column per allele and each consecutive column pair being the 2 alleles for one locus. Alleles are indicated with their microsatellite length. Adegenet files are either genepop (.gen) or arlequin files (.arp) and are used in the R package 'Adegenet' to do the majority of the genetic analyses performed in this study. The genepop files are similar to the format of data frames but use commas for population markers as well as a separate line with "POP" to indicate a population break and all individuals are named by their population ID. Also at the top of the document all the loci used to genotype the individuals are listed. There is also a comment line with what populations and individuals each adegenet file. The other 2 types of files included in the Data_Files folder are Structure_Files and Geneclass_Files. Structure_Files contains all the input and output files from the Structure runs performed in this study, and the Geneclass_Files contains all of the input and output files generated in the Geneclass2 analyses conducted for this study. 

<b>File Name Conventions</b>
<ul><li>Any file named as "clean" indicates that these files were cleaned in "2_clonecheck_md_relatedness.R" script, which cleans data frames and adegenet files for individuals with too much missing data (over 25%) and clones. These data files are used for all other genetic analyses in the study. </li></ul>
<ul><li>ZAIN_og files are the original data files used in analysis, with the originally scored genetic data. ZAIN_rebinned data files have been rebinned to make results more consistent by researcher (see explanation in "rebinning analysis" script). </li></ul>
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
<ul><ul><li>Garden_Wild_Comparisons</li></ul></ul>
<ul><ul><li>Scoring_Comparison</li></ul></ul>
<ul><ul><li>Sum_Stats</li></ul></ul>
