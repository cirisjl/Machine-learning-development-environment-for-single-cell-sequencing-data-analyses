# Refine Annotation for Mouse Bone Marrow using Visulization Tool[^1]

I this tutorial, we will show you how to use the Visulizaiton tool to refine the annotation mouse bone marrow single-cell sequencing data. Bone marrow was collected from uninfected WT and LCMV Cl 13-infected WT and Sphk2−/− mice using 10× Genomics single cell 3’ RNA-seq library [^1]. 

<hr/>

## Create a New Visulization Job

Go to **Analyses → Tools → Visulization** to start a new UMAP/t-SNE visulization job.

### Select Input Dataset(s)

Open the **Datasets** file selector. Use filters, global search and/or sorters to fast locate the tartget dataset **Mouse Marrow** with 19448 cell counts.

By clicking the **+** sign on the left of each column, you can expand the multi-stage data. We will choose the **CellTypist, scANVI, SingleR Annotation for SCVI Integration**.

<div align="center">
<img src="/images/doc/tool_visulization_input_dataset.png" width="80%">
</div>

> [!NOTE]  
> To use your own data, see **Doc → Data Management → Upload Data**.

### Choose Methods and Tune Parameters

Visualization parameters include `n_neighbors`, `n_pcs`, `Resolution` for UMAP/t-SNE and `Number of Highly Variable Genes for Heatmap`.

We will use **1** for `Resolution` and **50** for `Number of Highly Variable Genes for Heatmap`. The system will draw the top 50 marker genes for each cell type.

If `Cluster Column Name` is specified, it will be used for the annotation panel, and even UMAP/t-SNE is set to True, the new cluster will not be used for the annotation panel. Therefore, it's suggested to turn re-run UMAP/t-SNE off to save time. Otherwise, the default Leiden clustering will be used.

We will use **X_scVI_leiden** as `Cluster Column Name`. It's generated from scVI embeddings using Leiden clustering.

<div align="center">
<img src="/images/doc/tool_visulization_ccn" width="80%">
</div>

> [!TIP] 
> `n_neighbors` (Number of Neighbors): How many nearby cells are considered when building the graph for dimension reduction (e.g., UMAP) and clustering.
    + **Effect**: Small values emphasize **local details (more fragmented clusters/embedding)**; large values emphasize **global structure (more cohesive)**.
> `n_pcs` (Number of Principal Components): The number of dimensions kept after initial linear dimension reduction (PCA).
    + **Effect**: Too few can lose biological signal; too many can retain noise. Crucial for balancing signal and noise before further analysis.
> `Resolution`: A parameter for graph-based clustering (e.g., Leiden/Louvain) that controls cluster granularity.
    + **Effect**: Small values yield fewer, larger clusters; large values yield more, smaller clusters (more subtypes).

<hr/>

## Monitor Job

Once job is submitted, it will be queued and processed on backend. You can leave your current page and track your job status and retrieve the results on **"My Jobs"** page (see **"Job Management"**).

<div align="center">
<img src="/images/doc/tool_visulization_myjobs" width="80%">
</div>

The job details page will show you the **Dataset Information**, **Execution Details** and **Live Logs**, from which, you can see the excution details.

Ann error report panel will be shown if the job is failed. User can submit the issues to GitHub issues using this function. So, we can track and fix the issues.

<div align="center">
<img src="/images/doc/tool_visulization_job_details.png" width="80%">
</div>

> [!NOTE]  
> Single-cell research generates massive datasets. We optimize storage and processing by using a deduplication strategy. Each file is assigned a unique **MD5 hash** based on its **content, process, method, and parameters**. If a matching hash exists, we retrieve the precomputed result, saving significant time and computational resources.

## Results

All task outputs are provided in **`AnnData`** or **`MuData`** format. You can convert them into other formats anytime using the built-in **Formatting Tool (Toolos → Formatting → Convert)**.

### Visualization
Reduced-dimensional embeddings (e.g., UMAP, t-SNE, PCA) are stored in:

- `AnnData.obsm["layer_method"]`
- **Examples:**  
  - `AnnData.obsm["X_umap"]`  
  - `AnnData.obsm["sctransform_umap_3d"]`

User can access to interactive 2D/3D UMAP/t-SNE plots in the plots panel. You can change the **Color** dropdown list to change the labels of the plots.

<div align="center">
<img src="/images/doc/tool_visulization_process_umap.png" width="80%">
</div>

<div align="center">
<img src="/images/doc/tool_visulization_process_tsne.png" width="80%">
</div>

The **Gene Expression** panel shows **UMAP**, **Cell Sets**, **Cell Set Sizes**, **Heatmap** and **Gene List**. User can switch from different annotations to see the annotation, cell composition, and gene expression and use prior knowledge to determine if the annotation is reasonable.

<div align="center">
<img src="/images/doc/tool_visulization_process_ge.png" width="80%">
</div>

In the **Manual Annotaiton** panel, the system generate a table includes a majority vote (**cell_label**) from results of all methods, annotation of each method and the confident scores. Use can **Edit** and **Delete** the **Cluster**.

<div align="center">
<img src="/images/doc/tool_visulization_process_ma.png" width="80%">
</div>

> [!TIP]
> You may use the **Header** to sort and filter the data.
>
> If the required cell type is not in the dropdown list of **cell_type** column, you may type your cell type and create a new one.
>
> You may use the **Restore Deleted Rows** to restore the deleted rows.

## AI Assistant

The AI Assistant provide a way of seemless integrating prior knowledages from AI agents (e.g.m, chatGPT, Gemini) with the Tools and Workflows. For annotation, it passes the top 50 Highly Variable Genes (HVG) of each cluster to AI agents. User can also use the chatbox or upload paper to the chatbox to introduce more information to help with the annotation. User can download the chat history for review.

<div align="center">
<img src="/images/doc/tool_visulization_process_ai.png" width="80%">
</div>

<hr/>

## References
[^1]: Vijayamahantesh V, He Y, Jiang L, Huerter H, Jung KI, McKenna S, Studstill CJ, Allen L-AH, Nistala R, Xu D and Hahm B (2026) Sphingosine kinase 2 suppresses neutrophil responses to promote viral persistence while attenuating immune pathology. Front. Immunol. 16:1706967. doi: 10.3389/fimmu.2025.1706967