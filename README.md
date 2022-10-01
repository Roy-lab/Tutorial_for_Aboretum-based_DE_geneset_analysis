# Tutorial_for_Aboretum-based_DE_geneset_analysis
A tutorial which includes Arboretum &amp; findTransitionGenesets with an example dataset.

- [0. About this tutorial](#step-0-about-this-tutorial)
- [1. Prepare input tree](#step-1-prepare-input-tree)
- [2. Prepare input order](#step-2-prepare-input-order-file)
- [3. Prepare input OGID](#step-3-prepare-input-OGID-file)
- [4. Prepare input value](#step-4-prepare-input-value-files)
- [5. Run initialization clustering / prepare input config](#step-5-run-initialization-clustering-and-prepare-config-file)
- [6. Run Arboretum](#step-6-run-arboretum)
- [7. Run findTransitioningGenesets](#step-7-run-findTransitionGenesets)
- [8. Result visualization](#step-8-visualize-the-results)

---

### \[Step 0\] About this tutorial

This tutorial explains how to run the ***Arboretum*** and ***findTransitionGeneset***. Arboretum is a multi-task clustering framework which uses hierarchical relationship of samples simultaneously while grouping the genes into a finite number of expression states. In this tutorial, we will do the **clustering of genes** based on the **pseudo-bulk expression** of each dataset, i.e. average values of gene expressions within each dataset (cell cluster, in single cell setting). Once we have defined these **gene expression states** for all datasets, we can identify genes with interesting patterns of expression on the hierarchy and group the genes based on the similarity of their patterns to identify the **differential expressing (DE) genesets**. 

`input_files` folder contains tutorial dataset which consists of following files :

| Dataset | Description| File name | 
| :---    | :---  | :--- |
| c1_matrix | 155 cells x 20840 genes, with row & column headers | input_files/c1_matrix.txt |
| c2_matrix | 166 cells x 20840 genes, with row & column headers | input_files/c2_matrix.txt |
| c3_matrix | 126 cells x 20840 genes, with row & column headers | input_files/c3_matrix.txt |
| c4_matrix | 265 cells x 20840 genes, with row & column headers | input_files/c4_matrix.txt |
| c5_matrix | 155 cells x 20840 genes, with row & column headers | input_files/c5_matrix.txt |
| Gene ID | Gene names of 20840 genes | input_files/allgenenames.txt |
| Tree | [newick format](https://evolution.genetics.washington.edu/phylip/newicktree.html) tree text | input_files/newick_tree.txt 

- **Note**: the codes and scripts of this tutorial were fully tested and checked as working with the format above. So please firslty prepare your data as like above if you want to work with your own data.

`code` folder contains the main programs of this tutorial, which could be also downloaded from its own site:

| Program | Description | Original site |
| :---    | :---  | :--- |
| learnMoE | Initialization GMM clustering | https://github.com/Roy-lab/learnMoE |
| Arboretum | Multi-task gene clustering | https://github.com/Roy-lab/Arboretum2.0 |
| findTransitionGenesets | Hierarchical clustering of patterns identified from Arboretum | https://github.com/Roy-lab/clade-specific_gene_sets |
| reformatSpeciesTree | Format newick tree to Arboretum input | - |

- **Note**: All programs are written in C++ and the compling could be fulfilled by "Makefile" in each directory with the "make" command:


`script` folder contains several shell and PERL scripts which is used for the specific formatting of the data while running the wrapper scripts.

The key steps of clustering are performed by the following wrapper scripts:

| Wrapper script name | Running program | Output |
| :---    | :---  | :--- |
| step1_run_GMM_and_prep_config.sh | learnMoE | arb_input/merged_gmm_k#/, arb_input/config_k#.txt |
| step2_run_arboretum.sh | Arboretum | arb_output/ |
| step3_run_findTransitionGenesets.sh | findTransitionGenesets | transition_genesets_output/ |

---

### \[Step 1\] Prepare input tree

Arboretum requires a relationship tree structure between samples for the modeling. The format of the input tree is a text file consist of 3 columns (tab delimited). Each row is explaining the relationship of a parent node and a child node. One parental node always have 2 children nodes (left and right). The left/right children nodes will be dealt as equivalent, i.e. there is no order between children nodes. An ancestral node could be a child of another superordinate ancestral node. Following example is showing a tree consists of 3 nodes (c1, c2, c3) where c1 and c2 are more close. "Anc" means ancestral node.
```
c1 (TAB) left (TAB) Anc2
c2 (TAB) right (TAB) Anc2
Anc2 (TAB) left (TAB) Anc1
c3 (TAB) right (TAB) Anc1
```

To deal with a tree easier, especially for a complex tree with many leaves, we can use a program `reformatSpeciesTree` in the `code` folder. This program takes newick-style text file as a input and prints out the tree file for the Arboretum as the output. 

**Running command:**
```
- USAGE: reformatSpeciesTree [newick_tree.txt]

./code/reformatSpeciesTree input_files/newick_tree.txt > arb_input/tree.txt
```

**Result files:**
`input_files/newick_tree.txt` shows the newick tree of relationship among 5 sampples in this tutorial. The output `arb_input/tree.txt` is a formatted tree for the Arboretum.

---

### \[Step 2\] Prepare input order file

This text files is a simple list of sample names WITHOUT Ancestral nodes ("Anc#") with format like below:
```
c2
c4
c3
c1
c5
```

**Running command:**
```
cut -f1 arb_input/input_tree.txt |grep -v Anc > arb_input/orders.txt
```
This script used the tree file prepared in the previous step. 

**Result files:**
`arb_input/orders.txt` is the order file of our 5 samples example dataset, an input file of the Arboretum.

---

### \[Step 3\] Prepare input OGID file

OGID stands for "OrthoGroup ID" and it is a list of orthogroups(OGs) with a profiled list of corresponding gene IDs per species (wihch is a legacy of the original Arboretum). Regardless of the its naming, this file is used for the mapping of genes' identities from different datasets.
The format of this file is 
```
OGID (TAB) dataset1_GeneID,dataset2_GeneID,dataset3_GeneID,... 
```
(tab delimited between OGID and gene IDs, comma delimited among genes)

Note that the order of dataset-specific geneID should be same to the order file from [Step 2](#step-2-prepare-input-order-file). Also, from now on, the gene IDs of each dataset is named as "dataset_GeneID" (e.g. c1_AAA).

**Running command:**
```
- USAGE: generating_OGID.sh [orders.txt] [allgenenames.txt]

sh script/generating_OGID.sh arb_input/orders.txt input_files/allgenenames.txt > arb_input/OGID.txt
```
This script generates OGID based on the order of genes in "allgenenames.txt" file and shapes the gene IDs into "dataset#_geneID" as well as ordering it matched to the order file.

**Result files:**
The output `arb_input/OGID.txt` is an input file of the Arboretum.

---

### \[Step 4\] Prepare input value files

We are using pseudo-bulk expression vector of genes per samples, which is the average values of gene expressions within each dataset matrix so that each dataset has 1 vector values. The format of the file is a tab-delimited text file of 2 columns, the 1st column as "dataset_geneID" and the 2nd column is the mean expression values of all cells within the dataset. 
```
c1_A1BG	(TAB) 0.673579
c1_A1BG-AS1 (TAB) 0.246118
c1_A1CF (TAB) 0.000000
c1_A2M  (TAB) 0.899409
...
```

**Running command:**
```
- USAGE: generating_meanvals.pl [genes-by-cells_matrix.txt] ["dataset-name"]

perl script/generating_meanvals.pl input_files/c1_matrix.txt c1
mv c1_meanval.txt arb_input/
perl script/generating_meanvals.pl input_files/c2_matrix.txt c2
mv c2_meanval.txt arb_input/
... (Do this for all datset matrices.)
```

**Result files:** `arb_input/c#_meanval.txt` files are the values used as the source values for the Arboretum clustering.

**Transposing of the data matrix:** The script expects the matrix text file as a tab-delimited [genes x cells] matrix WITH row and column headers. If the user's dataset has given as [cells x genes], one can transpose it by using "script/transpose_matrix.pl" like:
```
- USAGE: transpose_matrix.pl [(tab-delimited) matrix.txt]
perl script/transpose_matrix.pl c1_matrix.txt > c1_matrix_transposed.txt
```

---

### \[Step 5\] Run initialization clustering and prepare config file

Identification of initializaiton clusters of genes is the first step of Arboretum analysis. This clustering is usually done by performing Gaussian Mixture Model (GMM) and this clustering assignments are assigned to the genes equally across all the datasets. Once after the initial clustering results are prepared, a speficific file for the configuration should be written for an input of Arboretum. The config file is a tab-delimited text file which has 3 columns:
```
DATASET_NAME (TAB) PATH/TO/INIT_CLUSTER.txt (TAB) PATH/TO/INPUT_VALUES.txt
```

Wrapper script `step1_run_GMM_and_prep_config.sh` is written to perform the initialization GMM smoothly as well as to turn the results as Arboretum input form automatically. Jobs done by this script are:
- Performing GMM clustering based on the merged matrix of input values in [Step 4](#step-4-prepare-input-value-files)
- Formatting the clustering result as dataset-sepcific initial cluster files
- Writing the config file automatically

**Running command:**
```
- USAGE: step1_run_GMM_and_prep_config.sh [#_cluster] [order.txt] [#_samples]

sh step1_run_GMM_and_prep_config.sh 3 arb_input/orders.txt 5
```

**Result files:** `arb_input/merged_gmm_k#/` contains the dataset-specific cluster initialization cluster assignment files. `arb_input/config.k#.txt` is the config file which is an input of Arboretum.

**About selection of the number of clusters (k):** The k number of GMM is usually matched to the k number of following Arboretum clustering, e.g. if we decided to use k=5 for the initialization, Arboretum should be also performed with k=5. However, if you specifically want to do the evaluation in GMM clusstering step, you can use the `likelihood.txt` in the GMM result folder (`arb_input/merged_gmm_k#/`).

---

### \[Step 6\] Run Arboretum

Arboretum identifies the expression states of genes which could be vary across different datasets. These states are driven by the expression amount of genes, which means higher number of cluster ID corresponds to the higher expression of the gene. **The increase of the number of states (i.e. number of clusters = k number) enables more precise patterning (e.g. 3 level differences in k=3 while 5-level differences in k=5), but also increases the complexity of the results as well.** 

Wrapper script `step2_run_arboretum.sh` is written to perform the Arboretum more easily by providing the proper input files. However, there are bunch of other parameters for the adjusted run of Arboretum. For the detailed usages of the parameters, please refer to the original [Arboretum github page](https://github.com/Roy-lab/Arboretum2.0). 

**Running command:**
```
- USAGE: step2_run_arboretum.sh [#_cluster] [order.txt] [OGID.txt] [tree.txt] [config.txt] ["representative_dataset_name"] [output_folder_name]

sh step2_run_arboretum.sh 3 arb_input/orders.txt arb_input/OGID.txt arb_input/input_tree.txt arb_input/config.txt c1 arb_output/
```
- **Note:** "representative_dataset_name" is for the specification of a specific dataset name for printing out the result text with the gene ID of that dataset. You can do not mind this much and just use any dataset name because or OGIDs cannot have "missing" gene IDs, unlike the multi-species analysis cases.

**Result files:** `arb_output/` folder contains all the results generated by the Arboretum run. Among the result file, `allspecies_clusterassign_lca_brk.txt` is the most summarized result table, which is a tab-delimited text file. The rows are gene IDs of the "representative_dataset_name" which was specified by user
and columns are all nodes (leaves and branches) in the tree. The table is sorted by the cluster assignment of Ancestor1 (root ancestor of the tree; positioned at the right-most column of the table). Dummy lines are the separators of different cluster assignments between two group of Ancestor1 clusters.

**About selection of the number of clusters (k):** The specification of k number of the Arboretum could be done in several ways, such as by using the calculaed "penalized likelihood" or "BIC/AIC" values in the `likelihood.txt` file in the result output folder `arb_output/`, or by calculating the "clustering evaluation metrics (e.g. SI, DBI, etc.)" using the clustering results and input values (as the merged matrix by horizontal concatenation of gene values of datasets.

**About random initialization:** The wrapper script runs Arboretum based on the random seeding setting (parameter swithch: "-r yes"), which means different runs of Arboretum will result in slightly differen results even with same usage of k number. Running several different times with same k setting and find the best result based on a certain criterion (e.g. clustering evaluation metric) could be helpful for overcoming local optima problem of the Arboretum's E-M algorithm.

---

### \[Step 7\] Run findTransitionGenesets

Arboretum stratifies the semi-continuous expression levels of genes of the datasets into finite numbers of expression states. This results in a gene's profile of expression states per datasets. The assessment based on these profiles is beneficial because it can catch the nuanced expressional changes which is not allowed in the DE analysis based on pairwise comparisons in usual.

We are interesting to identify the groups of genes that are sharing a common expression states profile because those gene sets may be working together as a pathway and they will generate a specific pattern across different datasets or biological contexts. To identify this type of genesets systematically, we do a hierarchical clustering based on the expression state profiles.

Wrapper script `step3_run_findTransitionGenesets.sh` is developed for performing the hierarchical clustering of profiles. For the detailed usage of the program `findTransitionGenesets`, please refer to the [clustering code github page](https://github.com/Roy-lab/clade-specific_gene_sets/tree/v2/de-novo_clustering). 

**Running command:**
```
- USAGE: step3_run_findTransitionGenesets.sh [Arboretum_output_folder] [order.txt] [OGID.txt] ["representative_dataset_name"] [output_folder_name]

sh step3_run_findTransitionGenesets.sh arb_output/ arb_input/orders.txt arb_input/OGID.txt c1 transition_genesets_output
```

Practically, the jobs which are done by this wrapper script are:
- Reordering the cluster IDs given by Arboretum into expression amount-related orders
- Re-format the result table of Arboretum `allspecies_clusterassign_lca_brk.txt` to the input form of the hierarchical clustering (`allcelltypes_clusterassign_brk.txt`)
- Run the hierarchical clustering with given input parameters

**Result files:**
- **all_assign.txt** : a formatted result file for plotting the all genes (rows) with expression state profiles and average expression values. The file is written as the input format for [Heatmap.awk](http://compbio.mit.edu/pouyak/software/Heatmap.awk).
- **ordered_clusterset_means.txt** : a [Heatmap.awk](http://compbio.mit.edu/pouyak/software/Heatmap.awk)-formatted result file for plotting the pattern of identified clustersets as rows.
- **clusterset##.txt** : a [Heatmap.awk](http://compbio.mit.edu/pouyak/software/Heatmap.awk)-formatted result files for plotting cluster-specific heatmaps with genes within the cluster as rows.
- **all_genes_clusterassignment_matrix.txt** : a text table of the hierarchical clustering results in plain text. Rows are genes and columns are nodes of the relationship tree. Each clusters are delimited by "Dummy" rows which are starting with "Dummy###".
- **ordered_mean_clusterassign_matrix.txt** : a text table of the hierarchical clustering results in plain text but the rows are the clusterset IDs.
- **all_genesets.txt** : a text file of the hierarchical clustering result in different format. The 2 columns of this files are "cluster ID" and "lists of genes" (delimited as #).

---

### \[Step 8\] Visualize the results

[Step 7](#step-7-run-findTransitionGenesets)'s following result files are specifically formatted for the heatmap plotting using the script [Heatmap.awk](http://compbio.mit.edu/pouyak/software/Heatmap.awk).
- `all_assign.txt`
- `ordered_clusterset_means.txt`
- `clusterset##.txt` 

Following script will draw the heatmap and save as `*.svg` files based on the default color scheme. Please check the svg files in `transition_genesets_output/` folder to check how they look.
- `all_assign.svg`
- `ordered_clusterset_means.svg`
- `clusterset##.svg` 

**Running command:**
```
- USAGE: draw_heatmap.sh [formatted_input_file] [max_value_of_cluster] [max_value_for_depicting_expreesions]

sh script/draw_heatmap.sh transition_genesets_output/ordered_clusterset_means.txt 3 2
sh script/draw_heatmap.sh transition_genesets_output/clusterset101.txt 3 2
```

