// Development
export const NODE_API_URL = `http://${process.env.REACT_APP_HOST_URL}:3001/node`;
export const CELERY_BACKEND_API = `http://${process.env.REACT_APP_HOST_URL}:5005/api`;
export const FLASK_BACKEND_API = `http://${process.env.REACT_APP_HOST_URL}:5003`;
export const WEB_SOCKET_URL = `ws://${process.env.REACT_APP_HOST_URL}:5005/wsapi`;
export const UPPY_API_URL = `http://${process.env.REACT_APP_HOST_URL}:3020`;
export const DIRECTUS_URL = `http://${process.env.REACT_APP_HOST_URL}:8055`;

// Production
// export const NODE_API_URL = `https://${process.env.REACT_APP_HOST_URL}/node`;
// export const CELERY_BACKEND_API = `https://${process.env.REACT_APP_HOST_URL}/api`;
// export const FLASK_BACKEND_API = `https://${process.env.REACT_APP_HOST_URL}:5003`;
// export const WEB_SOCKET_URL = `wss://${process.env.REACT_APP_HOST_URL}:5005/wsapi`;
// export const UPPY_API_URL = `https://${process.env.REACT_APP_HOST_URL}`;
// export const DIRECTUS_URL = `https://${process.env.REACT_APP_HOST_URL}`;

// Storage
export const PUBLIC_DATASETS = "/usr/src/app/storage/publicDatasets/";
export const STORAGE = "/usr/src/app/storage";
// Github
export const owner = process.env.REACT_APP_OWNER;
export const repo = process.env.REACT_APP_REPO;

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
  export const defaultQcParams = {
    assay: "RNA",
    min_genes: 200,
    max_genes: 20000,
    min_cells: 2,
    target_sum: 10000,
    n_top_genes: 2000,
    n_neighbors: 15,
    n_pcs: 1,
    resolution: 1,
    doublet_rate: 0,
    regress_cell_cycle: false
  };

  export const defaultNormalizationParams = {
    assay: "RNA",
    n_neighbors: 15,
    n_pcs: 1,
    resolution: 1,
  };

  export const defaultReductionParams = {
    assay: "RNA",
    n_neighbors: 15,
    n_pcs: 1,
    resolution: 1,
  };
