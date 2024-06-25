export const uiSchema = {
  "parameters": {
    "classNames": "category",
    "output_format": {
      "classNames": "sub-category",
      "ui:widget": "SelectComponent",
      'ui:options': {
        clearable: true ,
        placeholder: "Select the Output Format",
        creatable: false,
        searchable: true,
        opts:["AnnData", "SingleCellExperiment", "Seurat", "CSV"] 
      }
    },
    "methods": {
      "classNames": "sub-category",
      "ui:widget": "MultiSelectComponent",
    },
    "params": {
      "classNames": "form-subset sub-category",
      "default_assay": {
        "classNames": "sub-category",
        "ui:widget": "ClusterLabelInput"
      },
      "npcs": {
        "classNames": "sub-category",
        "ui:widget": "RangeSlider",
        'ui:title': 'n_pcs: ', 
        'ui:options': {
          title: 'n_pcs: ', 
          min: 0,
          max: 200,
          step: 1,
          marks:[
            { value: 0, label: '0*' },
            { value: 5, label: '5' },
            { value: 10, label: '10' },
            { value: 20, label: '20' },
            { value: 40, label: '40' },
            { value: 50, label: '50' },
            { value: 125, label: '125' },
            { value: 200, label: '200' },
          ]
        }
      },
      "layer": {
        "classNames": "sub-category",
        "ui:widget": "ClusterLabelInput"
      },
      "species": {
        "classNames": "sub-category",
        "ui:widget": "SelectComponent",
        'ui:options': {
          clearable: true ,
          placeholder: "Select the Species type",
          creatable: false,
          searchable: true,
          opts:["human", "mouse"]
        }
      },
      "idtype": {
        "classNames": "sub-category",
        "ui:widget": "SelectComponent",
        'ui:options': {
          clearable: true ,
          placeholder: "Select the ID type",
          creatable: false,
          searchable: true,
          opts:["SYMBOL", "ENSEMBL", "ENTREZID", "REFSEQ"]
        }
      },
      "min_genes": {
        "classNames": "sub-category",
        "ui:widget": "RangeSlider",
        'ui:options': {
          title: 'Min Genes: ', // Title for the slider
          min: 0,
          max: 20000,
          step: 25,
          marks: [
            { value: 200, label: '200*' },
            { value: 1000, label: '1000' },
            { value: 5000, label: '5000' },
            { value: 10000, label: '10000' },
            { value: 15000, label: '15000' },
            { value: 20000, label: '20000' }
          ]
        },
        'ui:title': 'Min Genes', 
      },
      "max_genes": {
        "classNames": "sub-category",
        "ui:widget": "RangeSlider",
        'ui:options': {
          title: 'Max Genes: ', // Title for the slider
          min: 0,
          max: 20000,
          step: 25,
          marks: [
            { value: 200, label: '200*' },
            { value: 1000, label: '1000' },
            { value: 5000, label: '5000' },
            { value: 10000, label: '10000' },
            { value: 15000, label: '15000' },
            { value: 20000, label: '20000' }
          ]
        },
        'ui:title': 'Max Genes', 
      },
      "show_umap": {
        "classNames": "sub-category",
        "ui:widget": "toggle"
      },
      "show_error": {
        "classNames": "sub-category",
        "ui:widget": "toggle"
      },
    }
  }
};