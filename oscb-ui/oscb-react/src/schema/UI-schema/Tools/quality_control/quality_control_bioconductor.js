export const uiSchema = (dynamicOptions) => ({
  "parameters": {
    "classNames": "category",
      // "output_format": {
      //   "classNames": "sub-category",
      //   "ui:widget": "SelectComponent",
      //   'ui:options': {
      //     clearable: true ,
      //     placeholder: "Select the Output Format",
      //     creatable: false,
      //     searchable: true,
      //     opts: dynamicOptions.outputFormatOptions 
      //   }
      // },
      "species": {
        "classNames": "sub-category",
        "ui:widget": "SelectComponent",
        'ui:options': {
          clearable: true ,
          placeholder: "Select the Species type",
          creatable: false,
          searchable: true,
          opts: dynamicOptions.speciesOptions
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
      "cluster_label": {
        "classNames": "sub-category",
        "ui:widget": "ClusterLabelInput"
      },
      // "show_umap": {
      //   "classNames": "sub-category",
      //   "ui:widget": "toggle"
      // },
      // "show_error": {
      //   "classNames": "sub-category",
      //   "ui:widget": "toggle"
      // },
      "qc_params": {
        "classNames": "form-subset sub-category",
        "methods": {
          "classNames": "sub-category",
          "ui:widget": "MultiSelectComponent",
        },
        "assay": {
          "classNames": "sub-category",
          "ui:widget": "ClusterLabelInput"
        },
        // "layer": {
        //   "classNames": "sub-category",
        //   "ui:widget": "SelectComponent",
        //   'ui:options': {
        //     clearable: true ,
        //     placeholder: "Select the Layer",
        //     creatable: false,
        //     searchable: true,
        //     opts: dynamicOptions.layers 
        //   }
        // },
        // "geneRange": {
        //   "classNames": "sub-category",
        //   "ui:widget": "geneRangeSlider",
        // },
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
        "min_cells": {
          "classNames": "sub-category",
          "ui:widget": "RangeSlider",
          'ui:options': {
            title: 'Min Cells: ', // Title for the slider
            min: 1,
            max: 200,
            step: 1,
            marks: [
              { value: 2, label: '2*' },
              { value: 10, label: '10' },
              { value: 50, label: '50' },
              { value: 100, label: '100' },
              { value: 200, label: '200' },
            ]
          },
          'ui:title': 'Min Cells', 
        },
        "target_sum": {
          "classNames": "sub-category",
          "ui:widget": "RangeSlider",
          'ui:title': 'Target Sum: ', 
          'ui:options': {
            title: 'Target Sum', // Title for the slider
            min: 0,
            max: 1e6,
            step: 1e4,
            marks: [
              { value: 1e4, label: '1e4*' },
              { value: 1e5, label: '1e5' },
              { value: 1e6, label: '1e6' },
            ]
          }
        },
        "n_top_genes": {
          "classNames": "sub-category",
          "ui:widget": "RangeSlider",
          'ui:title': 'Highly Variable Genes (n_top_genes):', 
          'ui:options': {
            title: 'Highly Variable Genes (n_top_genes)', // Title for the slider
            min: 100,
            max: 10000,
            step: 25,
            marks: [
              { value: 100, label: '100' },
              { value: 500, label: '500' },
              { value: 1000, label: '1000' },
              { value: 2000, label: '2000*' },
              { value: 5000, label: '5000' },
              { value: 10000, label: '10000'}]
          }
        },
        "n_neighbors": {
          "classNames": "sub-category",
          "ui:widget": "RangeSlider",
          'ui:title': 'n_neighbors: ', 
          'ui:options': {
            title: 'n_neighbors: ', // Title for the slider
            min: 2,
            max: 100,
            step: 1,
            marks:[
              { value: 2, label: '2' },
              { value: 5, label: '5' },
              { value: 10, label: '10' },
              { value: 15, label: '15*' },
              { value: 20, label: '20' },
              { value: 50, label: '50' },
              { value: 100, label: '100' },
            ]
          }
        },
        "n_pcs": {
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
        "resolution": {
          "classNames": "sub-category",
          "ui:widget": "RangeSlider",
          'ui:title': 'Resolution: ', 
          'ui:options': {
            title: 'Resolution: ', 
            min: 0,
            max: 5,
            step: 0.05,
            marks:[
              { value: 0.1, label: '0.1' },
              { value: 0.5, label: '0.5*' },
              { value: 1, label: '1' },
              { value: 2.5, label: '2.5' },
              { value: 5, label: '5' },
            ]
          }
        },
        "doublet_rate": {
          "classNames": "sub-category",
          "ui:widget": "RangeSlider",
          'ui:title': 'Expected Doublet Rate:', 
          'ui:options': {
            title: 'Expected Doublet Rate:', // Title for the slider
            min: 0,
            max: 0.5,
            step: 0.001,
            marks: [
              { value: 0, label: '0%' },
              { value: 0.08, label: '8%*' }, // Default
              { value: 0.125, label: '12.5%' },
              { value: 0.2, label: '20%' },
              { value: 0.5, label: '50%' },
            ]
          }
        },
        "regress_cell_cycle": {
          "classNames": "sub-category",
          "ui:widget": "SwitchComponent"
        },
        "colour_by": {
          "classNames": "sub-category",
          "ui:widget": "ClusterLabelInput"
        },
        "shape_by_1": {
          "classNames": "sub-category",
          "ui:widget": "ClusterLabelInput"
        },
        "shape_by_2": {
          "classNames": "sub-category",
          "ui:widget": "ClusterLabelInput"
        },
        "use_default": {
          "classNames": "sub-category",
          "ui:widget": "toggle"
        }
      }
  }
});
