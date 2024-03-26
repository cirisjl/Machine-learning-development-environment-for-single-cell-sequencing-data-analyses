export const LOGIN_API_URL = `http://${process.env.REACT_APP_HOST_URL}:3001`;
export const SERVER_URL = "http://" + process.env.REACT_APP_HOST_URL + ":3001";
export const PREVIEW_DATASETS_API = `http://${process.env.REACT_APP_HOST_URL}:3001`;
export const CELERY_BACKEND_API = `http://${process.env.REACT_APP_HOST_URL}:5000`;
export const FLASK_BACKEND_API = `http://${process.env.REACT_APP_HOST_URL}:5003`;
export const PUBLIC_DATASETS = "/usr/src/app/storage/publicDatasets/"
export const STORAGE = "/usr/src/app/storage";
export const defaultValues = {
    min_genes: 200,
    max_genes: 20000, // No limit
    min_cells: 2,
    target_sum: 1e4,
    n_top_genes: 2000,
    n_neighbors: 15,
    n_pcs: 0, // None
    resolution: 1,
    regress_cell_cycle: false,
    use_default: true,
    doublet_rate: 0.08
  };
