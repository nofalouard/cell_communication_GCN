# CellcommGCN

**Description:** A graph convolutional network for identifying closely interacting cell types and novel ligand and receptors in single cell data.

### Using CellcommGCN

1. Create and activate the conda environment with necessary packages:

    ```shell
    conda env create -f environment.yml
    conda activate gcn
    ```


2. Import single cell data:
Single cell datasets should be in '.h5ad' file format, i.e. Anndata files. In the cell_communication_GCN/sc_data folder, create a subfolder for each file or set of files you would like to analyze separately. For example, you could create subfolders that seperate your datasets based on tissue or disease of origin. 

    All resulting analyses for a particular subfolder will be placed in the output directory `results_{subfolder_name}` 


3. Run cellcommGCN for a selected subfolder:

    ```shell 
    python analyze_adj.py folder_name
    ```
    This will build a graph representation of each dataset from your specified folder and apply a graph convolutional layer on each to identify graph domains of strong cell interaction and coherent gene expression.

    Outputted figures for each dataset:

    - `cell_type_ratios/clust_composition_plt_{dataset_name}.pdf`: Barplot displaying the percentage of cells from each cell type lying in each one of a dataset's identified graph domains.

    - Intermediate results necessary to run steps 4 and 5 


4. Analyze interaction strength between cell types:

    ```shell 
    python cell_pair_hmap.py folder_name
    ```

    Using intermediate results from step 3, this will generate heatplots for each cluster identified in a dataset, displaying the interaction score between cells of every cell type pair. 

    The score is defined as the log transformed p value of the mean edge weight lying between cells of a given cell type pair, e.g. between macrophages and CD4 T cells. The significance is based on a null distribution of the mean edge weight, generated by randomly permuting all cell type labels in the dataset. Only significant interaction scores (p value < 0.05) are displayed on the heatmap.

    Outputted figures for each dataset:

    - `interaction_figures/{dataset_name}_{cluster_number}_Interactions.pdf` 



5. Identify significantly co-expressed genes among interacting cells 

    ```shell
    python analyze_coexpression.py folder_name
    ```

    This finds gene pairs that the cells of each graph domain (identified in step 3) co-express for each dataset stored under the specified folder. It then outputs the common gene pairs significantly co-expressed in all the datasets in a csv file. These co-expressed genes can potentially serve as ligand-receptor pairs used among closely interacting cells in their communication pathways.

    Outputted file:

    - `co_expression/coexprsd_genes_{foldername}.csv`











