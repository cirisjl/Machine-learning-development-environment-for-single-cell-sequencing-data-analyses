export const uiSchema = () => ({
  "parameters": {
    "ui:classNames": "category",
    "description": {
      "ui:classNames": "sub-category",
      "ui:widget": "ClusterLabelInput",
    },
    "mode": {
      "ui:classNames": "sub-category",
      "ui:widget": "SelectComponent",
      'ui:options': {
        clearable: false,
        creatable: false,
        searchable: true,
        opts: [
          "PHOENIX",
          "CDC_PHOENIX",
          "SCAFFOLD",
          "CDC_SCAFFOLD",
          "SRA",
          "CDC_SRA",
          // "combine_griphins"
        ]
      }
    },
    "kraken2db": {
      "ui:classNames": "sub-category",
      "ui:widget": "SelectComponent",
      'ui:options': {
        clearable: false,
        creatable: false,
        searchable: true,
        opts: [
          "Standard-8"
        ]
      }
    },
    "scaffolds_params": {
      "ui:classNames": "form-subset sub-category",
      "scaffold_ext": {
        "ui:classNames": "sub-category",
        "ui:widget": "ClusterLabelInput"
      },
    },
    "sra_params": {
      "ui:classNames": "form-subset sub-category",
      "use_sra": {
        "ui:classNames": "sub-category",
        "ui:widget": "toggle"
      },
    },
    "species_params": {
      "ui:classNames": "form-subset sub-category",
      "centar": {
        "ui:classNames": "sub-category",
        "ui:widget": "toggle"
      },
      "griphin_out": {
        "ui:classNames": "sub-category",
        "ui:widget": "ClusterLabelInput"
      },
    },
    "optional_params": {
      "ui:classNames": "form-subset sub-category",
      "coverage": {
        "ui:classNames": "sub-category",
        "ui:widget": "RangeSlider",
        // 'ui:title': 'If you want to increase the coverage cut off >30x use --coverage: ',
        'ui:options': {
          // title: 'If you want to increase the coverage cut off >30x use --coverage: ', 
          min: 1,
          max: 100,
          step: 1,
          marks: [
            { value: 20, label: '20' },
            { value: 30, label: '30*' },
            { value: 40, label: '40' },
            { value: 50, label: '50' },
            { value: 60, label: '60' },
            { value: 70, label: '70' },
            { value: 80, label: '80' },
            { value: 100, label: '100' },
          ]
        }
      },
      "create_ncbi_sheet": {
        "ui:classNames": "sub-category",
        "ui:widget": "toggle"
      },
      "email": {
        "ui:classNames": "sub-category",
        "ui:widget": "ClusterLabelInput"
      },
      "multiqc_title": {
        "ui:classNames": "sub-category",
        "ui:widget": "ClusterLabelInput"
      },
    },
  }
});