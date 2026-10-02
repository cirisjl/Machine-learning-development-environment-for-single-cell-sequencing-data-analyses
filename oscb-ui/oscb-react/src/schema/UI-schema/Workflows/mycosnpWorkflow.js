export const uiSchema = () => ({
  "parameters": {
    "ui:classNames": "category",
    "description": {
      "ui:classNames": "sub-category",
      "ui:widget": "ClusterLabelInput",
    },
    "workflow": {
      "ui:classNames": "sub-category",
      "ui:widget": "SelectComponent",
      'ui:options': {
        clearable: false,
        creatable: false,
        searchable: true,
        opts: [
          "NFCORE_MYCOSNP",
          "PRE_MYCOSNP",
        ]
      }
    },
    "pre_mycosnp_params": {
      "ui:classNames": "form-subset sub-category",
      "assembler": {
        "ui:classNames": "sub-category",
        "ui:widget": "SelectComponent",
        'ui:options': {
          clearable: false,
          creatable: false,
          searchable: true,
          opts: [
            "skesa",
            "velvet",
            "megahit",
            "spades"
          ]
        }
      },
      "shovill_depth": {
        "ui:classNames": "sub-category",
        "ui:widget": "RangeSlider",
        'ui:options': {
          min: 1,
          max: 100,
          step: 1,
          marks: [
            { value: 20, label: '20' },
            { value: 30, label: '30' },
            { value: 40, label: '40' },
            { value: 50, label: '50' },
            { value: 60, label: '60' },
            { value: 70, label: '70*' },
            { value: 80, label: '80' },
            { value: 100, label: '100' },
          ]
        }
      },
      "genome_size": {
        "ui:classNames": "sub-category",
        "ui:widget": "ClusterLabelInput"
      },
      "min_contig_cov": {
        "ui:classNames": "sub-category",
        "ui:widget": "RangeSlider",
        'ui:options': {
          min: 1,
          max: 100,
          step: 1,
          marks: [
            { value: 10, label: '10*' },
            { value: 20, label: '20' },
            { value: 30, label: '30' },
            { value: 40, label: '40' },
            { value: 50, label: '50' },
            { value: 60, label: '60' },
            { value: 80, label: '80' },
            { value: 100, label: '100' },
          ]
        }
      },
      "min_contig_len": {
        "ui:classNames": "sub-category",
        "ui:widget": "RangeSlider",
        'ui:options': {
          min: 100,
          max: 1000,
          step: 1,
          marks: [
            { value: 100, label: '100' },
            { value: 200, label: '200' },
            { value: 300, label: '300*' },
            { value: 400, label: '400' },
            { value: 500, label: '500' },
            { value: 600, label: '600' },
            { value: 800, label: '800' },
            { value: 1000, label: '1000' },
          ]
        }
      },
    },
    "nfcore_mycosnp_params": {
      "ui:classNames": "form-subset sub-category",
      "species": {
        "ui:classNames": "sub-category",
        "ui:widget": "SelectComponent",
        'ui:options': {
          clearable: false,
          creatable: true,
          searchable: true,
          opts: [
            "candida_auris_gca_016772135.1",
          ]
        }
      },
      "positions": {
        "ui:classNames": "sub-category",
        "ui:widget": "ClusterLabelInput"
      },
      "genes": {
        "ui:classNames": "sub-category",
        "ui:widget": "ClusterLabelInput"
      },
      "exclude": {
        "ui:classNames": "sub-category",
        "ui:widget": "ClusterLabelInput"
      },
      "save_reference": {
        "ui:classNames": "sub-category",
        "ui:widget": "toggle"
      },
      "save_alignment": {
        "ui:classNames": "sub-category",
        "ui:widget": "toggle"
      },
      "rapidnj": {
        "ui:classNames": "sub-category",
        "ui:widget": "toggle"
      },
      "fasttree": {
        "ui:classNames": "sub-category",
        "ui:widget": "toggle"
      },
      "snpeff": {
        "ui:classNames": "sub-category",
        "ui:widget": "toggle"
      },
      "iqtree": {
        "ui:classNames": "sub-category",
        "ui:widget": "toggle"
      },
      "raxmlng": {
        "ui:classNames": "sub-category",
        "ui:widget": "toggle"
      },
      "save_debug": {
        "ui:classNames": "sub-category",
        "ui:widget": "toggle"
      },
      "skip_samples": {
        "ui:classNames": "sub-category",
        "ui:widget": "ClusterLabelInput"
      },
      "skip_combined_analysis": {
        "ui:classNames": "sub-category",
        "ui:widget": "toggle"
      },
      "skip_phylogeny": {
        "ui:classNames": "sub-category",
        "ui:widget": "toggle"
      },
      "save_debug": {
        "ui:classNames": "sub-category",
        "ui:widget": "toggle"
      },
    },
    "sample_ploidy": {
      "ui:classNames": "sub-category",
      "ui:widget": "RangeSlider",
      'ui:options': {
        min: 1,
        max: 100,
        step: 1,
        marks: [
          { value: 1, label: '1*' },
          { value: 10, label: '10' },
          { value: 20, label: '20' },
          { value: 30, label: '30' },
          { value: 40, label: '40' },
          { value: 50, label: '50' },
          { value: 60, label: '60' },
          { value: 100, label: '100' },
        ]
      }
    },
    "coverage": {
      "ui:classNames": "sub-category",
      "ui:widget": "RangeSlider",
      'ui:options': {
        min: 0,
        max: 100,
        step: 1,
        marks: [
          { value: 0, label: '0*' },
          { value: 20, label: '20' },
          { value: 30, label: '30' },
          { value: 40, label: '40' },
          { value: 50, label: '50' },
          { value: 60, label: '60' },
          { value: 70, label: '70' },
          { value: 100, label: '100' },
        ]
      }
    },
    "gvcfs_filter": {
      "ui:classNames": "sub-category",
      "ui:widget": "ClusterLabelInput"
    },
    "gatkgenotypes_filter": {
      "ui:classNames": "sub-category",
      "ui:widget": "ClusterLabelInput"
    },
    "max_amb_samples": {
      "ui:classNames": "sub-category",
      "ui:widget": "RangeSlider",
      'ui:options': {
        min: 1000000,
        max: 100000000,
        step: 1000000,
        marks: [
          { value: 10000000, label: '10M*' },
          { value: 20000000, label: '20M' },
          { value: 30000000, label: '30M' },
          { value: 40000000, label: '40M' },
          { value: 50000000, label: '50M' },
          { value: 60000000, label: '60M' },
          { value: 80000000, label: '80M' },
          { value: 100000000, label: '100M' },
        ]
      }
    },
    "max_perc_amb_samples": {
      "ui:classNames": "sub-category",
      "ui:widget": "RangeSlider",
      'ui:options': {
        min: 0,
        max: 100,
        step: 1,
        marks: [
          { value: 0, label: '0' },
          { value: 10, label: '10*' },
          { value: 20, label: '20' },
          { value: 30, label: '30' },
          { value: 40, label: '40' },
          { value: 50, label: '50' },
          { value: 60, label: '60' },
          { value: 80, label: '80' },
          { value: 100, label: '100' },
        ]
      }
    },
    "min_depth": {
      "ui:classNames": "sub-category",
      "ui:widget": "RangeSlider",
      'ui:options': {
        min: 0,
        max: 100,
        step: 1,
        marks: [
          { value: 0, label: '0' },
          { value: 10, label: '10*' },
          { value: 20, label: '20' },
          { value: 30, label: '30' },
          { value: 40, label: '40' },
          { value: 50, label: '50' },
          { value: 60, label: '60' },
          { value: 80, label: '80' },
          { value: 100, label: '100' },
        ]
      }
    },
    "mask": {
      "ui:classNames": "sub-category",
      "ui:widget": "toggle"
    },
  }
});