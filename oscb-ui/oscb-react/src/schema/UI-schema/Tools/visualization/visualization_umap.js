export const uiSchema = (dynamicOptions) => ({

    "parameters": {
      "classNames": "category",
        "cluster_label": {
          "classNames": "sub-category",
          "ui:widget": "ClusterLabelInput"
        },
        "reduction_params": {
          "classNames": "form-subset sub-category",
          "assay": {
            "classNames": "sub-category",
            "ui:widget": "ClusterLabelInput"
          },
          "use_rep": {
            "classNames": "sub-category",
            "ui:widget": "SelectComponent",
            'ui:options': {
              clearable: false,
              placeholder: "Use the indicated representation",
              creatable: false,
              searchable: true,
              opts: dynamicOptions.embeddings
            }
          },
          "layer": {
            "classNames": "sub-category",
            "ui:widget": "SelectComponent",
            'ui:options': {
            clearable: false,
            placeholder: "Select the Layer",
            creatable: false,
            searchable: true,
            opts: dynamicOptions.layers 
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
                { value: 0, label: '0' },
                { value: 5, label: '5' },
                { value: 10, label: '10' },
                { value: 20, label: '20*' },
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
          "use_default": {
            "classNames": "sub-category",
            "ui:widget": "toggle"
          }
        }
    }
  });