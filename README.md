# Tutorial_for_Aboretum-based_DE_geneset_analysis
A tutorial which includes Arboretum &amp; findTransitionGenesets with an example dataset.

- [0. About this tutorial](#step-0-about-this-tutorial)
- [1. Install the program](#step-1-install-the-program)
- [2. Prepare input files](#step-2-prepare-input-data)
- [3. Prepare input tree](#step-3-prepare-input-tree)
- [4. Run TMF](#step-4-run-tmf)
- [5. Generate output clusters](#step-5-generate-output-clusters)

### \[Step 0\] About this tutorial

This tutorial explains how to run the ***Arboretum*** and ***findTransitionGeneset*** based on the result from Arboretum. Arboretum was originally developed for the gene exprssion analysis of multiple species, but the multi-task framework could be also used in other types of comparative analysis including comparative network anslsys or analysis of different cell groups.<br>
In this tutorial, we are assuming that we wanted to identify differential expressing geneset which exhibits complex patterns across multiple samples. Before going into the procedure, Please check the following initial files are correctly existing.

-  [ ] `input_files` folder : contains 5 different data matrices where the cells are different but the genes are common acorss all samples. This folder also contains newick tree file (newick.txt) for the relationship tree 


### \[Step 1\] Install the program

Follow the instruction at the [TMF github page](https://github.com/Roy-lab/tmf).

To check if the program has been installed properly, run the following commands (in shell) and see if the usage is printed out correctly without error messages.
```
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:PATH/TO/GSL/gsl-2.6/lib

./run_tmf -h
```


### \[Step 2\] Prepare input data

Check the following input files are prepared exactly as described.
- [ ] `samplenames.txt`: A text file of listing sample IDs (e.g. sampleX, sampleY, sampleZ, etc.) in a column.
- [ ] `allgenenames.txt` : A text file of all gene names in a column.
- [ ] `sampleX_cellnames.txt`,`sampleY_cellnames.txt`, ... : Text files of all cell names (barcodes) per samples.
- [ ] `sampleX_datamat.txt`,`sampleY_datamat.txt`, ... : TAB-delimited text data value matrix files per samples **WITHOUT** row and column header lines.

### \[Step 3\] Prepare input tree


### \[Step 4\] Run TMF


### \[Step 5\] Generate output clusters


See [Parameters](#parameters) for more details.
```
./run_tmf input/tree/toy_tree.txt 120 2 -o output/ -a 10 -l 200
```
- `input/toy_tree.txt` specifies the tree file, which contains file locations to individual task matrices (paths are relative to location of run_tmf executable location). 
- `120` is the number of features/columns in each task matrix, which has to be be the same across all tasks. 
- `2` = k, the smaller dimensions of U and V. 
-	[Optional] `-o output/` will put all output files to output/ directory. Check out the example output directory in the repo. By default output will be saved to current directory.
-	[Optional] `-a 10` will set the alpha (strength of regularization to parent node) to be 10. Default is alpha = 10.
- [Optional] `-l 200` will set lambda (strength of sparsity constraint) to be 200. By default there is no sparsity constraint, i.e., lambda = 0.

#### Input tree file format
See example in input/toy_tree.txt.
```
1 3 A input/toy/A.txt 95
2 3 B input/toy/B.txt 80
3 -1 root N/A N/A
```
- Column 1: **node ID**; start from 1 and move up.
- Column 2: **parent node ID**; right now the implementation will only work correctly if you ID all children nodes before a parent node (so start from the lowest depth of tree, move to next level, till you hit the root, which should be the last node ID.)
- Column 3: **node alias**, used as prefix to U and V output files.
- Column 4: **location of input matrix file for leaf nodes**, relative to where the `run_tmf` executable is. Set as N/A for non-leaf nodes.
- Column 5: **number of rows/data points in each input matrix**. Set as N/A for non-leaf nodes.
