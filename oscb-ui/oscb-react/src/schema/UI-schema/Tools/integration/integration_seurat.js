export const uiSchema = {
  "parameters": {
    "classNames": "category",
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
          marks: [
            { value: 0, label: '0' },
            { value: 5, label: '5' },
            { value: 10, label: '10' },
            { value: 20, label: '20' },
            { value: 30, label: '30*' },
            { value: 50, label: '50' },
            { value: 125, label: '125' },
            { value: 200, label: '200' },
          ]
        }
      },
      "dims": {
        "classNames": "sub-category",
        "ui:widget": "RangeSlider",
        'ui:title': 'dims: ',
        'ui:options': {
          title: 'dims: ', // Title for the slider
          min: 0,
          max: 200,
          step: 1,
          marks: [
            { value: 0, label: '0' },
            { value: 5, label: '5' },
            { value: 10, label: '10' },
            { value: 20, label: '20' },
            { value: 30, label: '30*' },
            { value: 50, label: '50' },
            { value: 125, label: '125' },
            { value: 200, label: '200' },
          ]
        },
      },
    }
  }
};