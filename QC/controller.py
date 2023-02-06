from fastapi import FastAPI
import uvicorn
import scanpy as sc

from services import processScanpyQC
from services import processDropkickQC

app = FastAPI()
    
@app.post('/processScanpyQC')
def processScanpy(filePath: str, min_genes: int, min_cells: int, n_genes_by_counts:int, pct_counts_mt:int):
    task = processScanpyQC.delay(filePath, min_genes, min_cells, n_genes_by_counts, pct_counts_mt)
    return {'taskId': task.task_id}

@app.post('/processDropkick')
def processDropkick(filePath: str, min_genes: int, min_cells: int, n_genes_by_counts:int, pct_counts_mt:int):
    task = processDropkickQC.delay(filePath, min_genes, min_cells, n_genes_by_counts, pct_counts_mt)
    return {'taskId': task.task_id}

if __name__ == '__main__':
    print(getattr(sc.datasets, 'pbmc68k_reduced'))
    uvicorn.run(app)
