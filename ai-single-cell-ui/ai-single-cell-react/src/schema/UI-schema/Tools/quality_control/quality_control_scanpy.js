export const uiSchema = {
    "parameters": {
      "classNames": "category",
        "output_format": {
          "classNames": "sub-category",
          "ui:widget": "select",
          "ui:placeholder": "Select file format"
        },
        "methods": {
          "classNames": "sub-category",
          "ui:placeholder": "Select a method",
          'ui:widget': () => (
            <div className='common-row-wrap'>
              <select>
                <option value="Scanpy">Scanpy</option>
              </select>
        </div>
          ),
        },
        "qualityControl": {
          "classNames": "category",
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
          "use_default": {
            "classNames": "sub-category",
            "ui:widget": "SwitchComponent"
          }
        },
        "projectionParameters": {
          "classNames": "category",
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
          }
        },
        "clusteringParameters": {
          "classNames": "category",
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
                { value: 0.5, label: '0.5' },
                { value: 1, label: '1*' },
                { value: 2.5, label: '2.5' },
                { value: 5, label: '5' },
              ]
            }
          }
        },
        "default_assay": {
          "classNames": "sub-category",
          'ui:widget': () => (
            <div className='common-row-wrap'>
              <span data-v-22825496="" class="ui-form-title-message warning"> * Optional </span>
              <input type='text' />
        </div>
          ),
        },
        "layer": {
          "classNames": "sub-category"
        },
        "path_of_scrublet_calls": {
          "classNames": "sub-category"
        },
        "species": {
          "classNames": "sub-category",
          "ui:placeholder": "Select species type"
        },
        "idtype": {
          "classNames": "sub-category"
        },
        "genes": {
          "classNames": "sub-category",
        },
        "ncores": {
          "classNames": "sub-category",
          "ui:widget": "range",
        },
        "show_umap": {
          "classNames": "sub-category",
          "ui:widget": "toggle"
        },
        "show_error": {
          "classNames": "sub-category",
          "ui:widget": "toggle"
        }
    }
  };
  