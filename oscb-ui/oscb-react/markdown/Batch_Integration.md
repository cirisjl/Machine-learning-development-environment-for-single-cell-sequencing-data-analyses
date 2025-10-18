## Task info
<div align="center">
<img src="http://c4130-110133.wisc.cloudlab.us:3000/images/evaluator/Batch Integration.png">
</div>

Single-cell RNA sequencing (scRNA-seq) experiments often face a pervasive challenge known as the "batch effect," where technical variations introduced during sample collection, preparation, or sequencing can create artificial differences between datasets (batches), even when the underlying biological samples are similar. As depicted in the figure's "Unintegrated Data" (Panel 1), cells of the same biological type (e.g., T Cells, represented by circles) from different batches (e.g., Batch 1 in blue, Batch 2 in red) appear as distinct, non-overlapping groups, while different cell types from the same batch might erroneously appear closer than they truly are. This technical noise complicates direct comparisons and downstream analyses, hindering the identification of true biological cell types and states. Single-cell batch integration methods, therefore, are critical computational techniques (illustrated as "Integration Magic" in Panel 2) designed to harmonize these diverse datasets. These methods work by identifying and aligning shared biological signals across different batches, effectively removing technical biases while preserving genuine biological distinctions. The outcome is "Integrated Data" (Panel 3), where cells of the same biological type (e.g., all T Cells, regardless of original batch) are correctly grouped together, allowing for accurate identification of cell populations and robust comparative analyses across multiple experiments.

<hr/>

## Metrics

This task evaluates batch integration methods by assessing their success in removing technical batch effects while preserving true biological variation (including ARI, ASW label, Cell Cycle Conservation, cLISI, HVG overlap, Isolated label F1 score, Isolated label ASW, and NMI) within single-cell data. Methods are provided with multi-batch, consistently labeled data (either normalized or unnormalized) and produce either a feature matrix, low-dimensional embedding, or neighborhood graph as integrated output. This output is subsequently evaluated using specific metrics that quantify both batch effect removal and biological variance conservation, with the task framework drawing from the most recent and comprehensive single-cell data integration benchmark.

A workflow for creating imputation benchmarks is available on [single-cell.ai](https://www.single-cell.ai/), where training data is stored in `adata.obsm["train"]` and test data in `adata.obsm["test"]`.

* **ARI:** The Adjusted Rand Index [^1] quantifies the congruence between a derived clustering solution and a predefined ground truth (e.g., cell types), factoring in chance agreement. It delivers a score between 0 (random assignment) and 1 (perfect concordance), reflecting both accurate inclusions and exclusions across partitions.

* **ASW batch:** This metric assesses the Poisson log-likelihood [^1] of observing the true gene counts in the test dataset, given that the denoised counts from the model represent the expected mean parameters of the underlying Poisson distributions.

* **ASW Label:** This metric assesses the Poisson log-likelihood [^1] of observing the true gene counts in the test dataset, given that the denoised counts from the model represent the expected mean parameters of the underlying Poisson distributions.

* **Cell Cycle Conservation:** This metric assesses the Poisson log-likelihood [^1] of observing the true gene counts in the test dataset, given that the denoised counts from the model represent the expected mean parameters of the underlying Poisson distributions.

* **cLISI:** This metric assesses the Poisson log-likelihood [^1] of observing the true gene counts in the test dataset, given that the denoised counts from the model represent the expected mean parameters of the underlying Poisson distributions.

* **Graph Connectivity:** This metric assesses the Poisson log-likelihood [^1] of observing the true gene counts in the test dataset, given that the denoised counts from the model represent the expected mean parameters of the underlying Poisson distributions.

* **HVG overlap:** This metric assesses the Poisson log-likelihood [^1] of observing the true gene counts in the test dataset, given that the denoised counts from the model represent the expected mean parameters of the underlying Poisson distributions.

* **iLISI:** This metric assesses the Poisson log-likelihood [^1] of observing the true gene counts in the test dataset, given that the denoised counts from the model represent the expected mean parameters of the underlying Poisson distributions.

* **Isolated label ASW:** This metric assesses the Poisson log-likelihood [^1] of observing the true gene counts in the test dataset, given that the denoised counts from the model represent the expected mean parameters of the underlying Poisson distributions.

* **Isolated label F1 score:** This metric assesses the Poisson log-likelihood [^1] of observing the true gene counts in the test dataset, given that the denoised counts from the model represent the expected mean parameters of the underlying Poisson distributions.

* **kBET:** This metric assesses the Poisson log-likelihood [^1] of observing the true gene counts in the test dataset, given that the denoised counts from the model represent the expected mean parameters of the underlying Poisson distributions.

* **NMI:** Normalized Mutual Information [^4] gauges the statistical dependence between a predicted clustering and a known categorical variable, normalizing by the average entropy of both distributions. This metric ranges from 0 (statistical independence) to 1 (perfect correlation), indicating the shared information content between the two partitions.

* **PCR:** Normalized Mutual Information [^4] gauges the statistical dependence between a predicted clustering and a known categorical variable, normalizing by the average entropy of both distributions. This metric ranges from 0 (statistical independence) to 1 (perfect correlation), indicating the shared information content between the two partitions.

<hr/>

## Python API

**Open Single-Cell Benchmarks (OSCB)** is a Python package with a collection of benchmark datasets, data loaders, and evaluators for single-cell machine learning.
To install OSCB, please use Python's package manager pip:

```bash
pip install oscb
```

To show the installed version, run:

```bash
python -c "import oscb; print(oscb.__version__)" 
```

You may update the version by running:

```bash
pip install -U oscb
```

OSCB's primary advantages are its intuitive **data loaders**, designed for seamless dataset ingestion, and its standardized **evaluators**, which enable consistent and objective performance assessment across imputation methods.

### Data loaders
To load a dataset, please replace `Benchmarks_ID/Dataset_ID` with the Benchmarks ID or Dataset ID (e.g., `"IM-h-10x_PBMC-1087-10x-2017"`). The default download folder is `./datasets`. You can change the data folder by passing a new value to the parameter `data_folder`.

For public Benchmarks and datasets, the downloaded file is in [AnnData](https://github.com/scverse/anndata) or [MuData](https://github.com/scverse/mudata) format. The raw counts are kept in `adata.X`. The processed data (normalized, imputed, ...) are stored in `adata.layers`.
```python
from oscb.data import DataLoader

# Download and process data at './datasets/'
adata = DataLoader("Benchmarks_ID/Dataset_ID", data_folder="./datasets/")
```

### Evaluators
#### OSCB Benchmarks
For Benchmarks procided by [single-cell.ai](https://www.single-cell.ai/), please replace `Benchmarks_ID` with the Benchmarks ID (e.g., `"IM-h-10x_PBMC-1087-10x-2017"`), and `denoised` with your denoised data (e.g., `adata.layer['denosied']` or `adata.obsm['denosied']`). To modify the method name, adjust the `method` parameter.
> [!IMPORTANT]
> The input types of `benchmarks_id` and `method` are `string`, while the input type of `denoised` is `{array-like, sparse matrix} of shape (n_samples, n_features)`.
```python
from oscb.evaluator import eval, write_json

results_dict = eval(adata, benchmarks_id="Benchmarks_ID", denoised=adata.layer['denosied'], method="Your method")
```
> [!TIP]  
> If a `benchmarks_id` is specified, OSCB will automatically generate a bar chart to visually compare the performance of the user's method against established benchmark approaches.
> <div align="center">
> <img src="http://c4130-110133.wisc.cloudlab.us:3000/images/evaluator/imputation_evaluation.png">
> </div>

#### User's datasets
To utilize `eval()`, whether loading a dataset within [single-cell.ai](https://www.single-cell.ai/) via its ID or using your own, you must supply the `task` type (e.g., `"Imputation"`  or `"IM"` for short), and an `denoised` (e.g., `adata.layer['denosied']` or `adata.obsm['denosied']`); the method name can be customized via the `method` parameter.
> [!IMPORTANT]
> The input types of `task` is `string`, while the input type of `denoised` is `{array-like, sparse matrix} of shape (n_samples, n_features)`.
```python
from oscb.evaluator import eval, write_json

results_dict = eval(task='Imputation', denoised=adata.layer['denosied'], method="Your method")
```
#### Save results
To save your results to JSON format file, please run:
```python
from oscb.evaluator import write_json

write_json(results, file_path="./output.json") # The default file path is ./output.json
```
#### Computing assessment
OSCB further offers a computational assessment unit that, by encapsulating your code within a `monitor` instance, tracks CPU, memory, GPU, and GPU memory usage.
```python
from oscb.utilization import Monitor

monitor = Monitor(1)
# Your code
monitor.stop()
```
<div align="center">
<img src="http://c4130-110133.wisc.cloudlab.us:3000/images/evaluator/imputation_utilization.png">
</div>

<hr/>

## References
[^1]: Hubert, L., & Arabie, P. (1985). Comparing partitions. Journal of Classification, 2(1), 193–218. https://doi.org/10.1007/bf01908075

[^4]: Luecken, M. D., Büttner, M., Chaichoompu, K., Danese, A., Interlandi, M., Mueller, M. F., Strobl, D. C., Zappia, L., Dugas, M., Colomé-Tatché, M., & Theis, F. J. (2021). Benchmarking atlas-level data integration in single-cell genomics. Nature Methods, 19(1), 41–50. https://doi.org/10.1038/s41592-021-01336-8