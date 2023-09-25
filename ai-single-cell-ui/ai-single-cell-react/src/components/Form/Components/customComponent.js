import React, { Component } from 'react';
import CreatableSelect from 'react-select/creatable';
import { SERVER_URL } from '../../../constants/declarations';

class MyForm extends Component {
  constructor(props) {
    super(props);
    this.state = {
      formData: {
        Dataset: '',
        Task: null,
        Downloads: '',
      },
      errors: {},
      isLoading: false,
      options: {
        Task: [], 
        Author: [],
        Species: [],
        'Sample Type':[],
        'Anatomical Entity': [],
        'Organ Part': [],
        'Model Organ': [],
        'Selected Cell Types': [],
        'Library Construction Method': [],
        'Nucleic Acid Source': [],
        'Disease Status (Specimen)': [],
        'Disease Status (Donor)': [],
        'Development Stage': [],
        'Cell Count Estimate': [],
        'Source': [],
      },
    };
  }

  componentDidMount() {
    // Make an API call to get the default options for all fields
    this.fetchDefaultOptions();
  }

  async fetchDefaultOptionsxx() {
    try {
      const response = await fetch(`${SERVER_URL}/mongoDB/api/options`);
      if (response.ok) {
        const data = await response.json();

        this.setState({
          options: {
            Task: data.Task.map((option) => ({ value: option, label: option })),
            Author: data.Author.map((option) => ({ value: option, label: option })),
            Species: data.Species.map((option) => ({ value: option, label: option })),
            'Sample Type': data['Sample Type'].map((option) => ({ value: option, label: option })),
            'Anatomical Entity': data['Anatomical Entity'].map((option) => ({ value: option, label: option })),
            'Organ Part': data['Organ Part'].map((option) => ({ value: option, label: option })),
            'Model Organ': data['Model Organ'].map((option) => ({ value: option, label: option })),
            'Selected Cell Types': data['Selected Cell Types'].map((option) => ({ value: option, label: option })),
            'Library Construction Method': data['Library Construction Method'].map((option) => ({ value: option, label: option })),
            'Nucleic Acid Source': data['Nucleic Acid Source'].map((option) => ({ value: option, label: option })),
            'Disease Status (Specimen)': data['Disease Status (Specimen)'].map((option) => ({ value: option, label: option })),
            'Disease Status (Donor)': data['Disease Status (Donor)'].map((option) => ({ value: option, label: option })),
            'Development Stage': data['Development Stage'].map((option) => ({ value: option, label: option })),
            'Cell Count Estimate': data['Cell Count Estimate'].map((option) => ({ value: option, label: option })),
            'Source': data['Source'].map((option) => ({ value: option, label: option })),
          },
        });
      } else {
        console.error('Error fetching default options');
      }
    } catch (error) {
      console.error('Error fetching default options:', error);
    }
  }

  async fetchDefaultOptions() {
    try {
      const response = await fetch(`${SERVER_URL}/mongoDB/api/options`);
      if (response.ok) {
        const data = await response.json();
  
        const options = {};
  
        // Check if each field exists in the response before mapping it
        if (data.Task) {
          options.Task = data.Task.map((option) => ({ value: option, label: option }));
        }
        if (data.Author) {
            options.Author = data.Author.map((option) => ({ value: option, label: option }));
        }
  
        if (data.Species) {
          options.Species = data.Species.map((option) => ({ value: option, label: option }));
        }
  
        if (data['Sample Type']) {
          options['Sample Type'] = data['Sample Type'].map((option) => ({ value: option, label: option }));
        }
  
        if (data['Anatomical Entity']) {
          options['Anatomical Entity'] = data['Anatomical Entity'].map((option) => ({ value: option, label: option }));
        }
  
        if (data['Organ Part']) {
          options['Organ Part'] = data['Organ Part'].map((option) => ({ value: option, label: option }));
        }
        if (data['Model Organ']) {
            options['Model Organ'] = data['Model Organ'].map((option) => ({ value: option, label: option }));
        }
        if (data['Selected Cell Types']) {
            options['Selected Cell Types'] = data['Selected Cell Types'].map((option) => ({ value: option, label: option }));
        }
        if (data['Library Construction Method']) {
            options['Library Construction Method'] = data['Library Construction Method'].map((option) => ({ value: option, label: option }));
        }
        if (data['Nucleic Acid Source']) {
            options['Nucleic Acid Source'] = data['Nucleic Acid Source'].map((option) => ({ value: option, label: option }));
        }
        if (data['Disease Status (Specimen)']) {
            options['Disease Status (Specimen)'] = data['Disease Status (Specimen)'].map((option) => ({ value: option, label: option }));
        }
        if (data['Disease Status (Donor)']) {
            options['Disease Status (Donor)'] = data['Disease Status (Donor)'].map((option) => ({ value: option, label: option }));
        }
        if (data['Development Stage']) {
            options['Development Stage'] = data['Development Stage'].map((option) => ({ value: option, label: option }));
        }
        if (data['Cell Count Estimate']) {
            options['Cell Count Estimate'] = data['Cell Count Estimate'].map((option) => ({ value: option, label: option }));
        }
        if (data['Source']) {
            options['Source'] = data['Source'].map((option) => ({ value: option, label: option }));
        }
  
        // Add other fields here
  
        this.setState({
          options,
        });
      } else {
        console.error('Error fetching default options');
      }
    } catch (error) {
      console.error('Error fetching default options:', error);
    }
  }
  

  handleChange = (e) => {
    const { name, value } = e.target;
    this.setState((prevState) => ({
      formData: {
        ...prevState.formData,
        [name]: value,
      },
    }));
  };

  handleTaskChange = (selectedOption) => {
    this.handleSelectChange('Task', selectedOption);
  };

  handleSelectChange(fieldName, selectedOption) {
    this.setState((prevState) => ({
      formData: {
        ...prevState.formData,
        [fieldName]: selectedOption,
      },
    }));
  }

  handleCreateOption = (fieldName, inputValue) => {
    this.setState((prevState) => ({
        formData: {
          ...prevState.formData,
          [fieldName]: inputValue,
        },
    }));
    this.setState((prevState) => {
      const newOption = { value: inputValue, label: inputValue };
      const updatedOptions = { ...prevState.options };
      updatedOptions[fieldName] = [...(updatedOptions[fieldName] || []), newOption];
      return {
        options: updatedOptions,
      };
    });
  };

  handleSubmit = (e) => {
    e.preventDefault();
    const errors = this.validateForm(this.state.formData);
    this.setState({ errors });

    if (Object.keys(errors).length === 0) {
      const formData = this.state.formData;
      console.log(formData);
    }
  };

  validateForm(formData) {
    const errors = {};

    if (!formData.Dataset) {
      errors.Dataset = 'Dataset is required';
    }

    if (!formData.Task || (formData.Task && formData.Task.value === '')) {
      errors.Task = 'Task is required';
    }

    if (!formData.Downloads) {
      errors.Downloads = 'Downloads is required';
    }

    return errors;
  }

  render() {
    const { formData, errors, isLoading, options } = this.state;
    return (
      <div>
        <h2>My Form</h2>
        <form onSubmit={this.handleSubmit}>
          {/* Dataset */}
          <div>
            <label>Dataset:</label>
            <input
              type="text"
              name="Dataset"
              value={formData.Dataset}
              onChange={this.handleChange}
            />
            {errors.Dataset && <p className="error">{errors.Dataset}</p>}
          </div>

          {/* Task (CreatableSelect) */}
          <div>
            <label>Task:</label>
            <CreatableSelect
              name="Task"
              value={formData.Task}
              isClearable
              isSearchable
              isLoading={isLoading}
              onChange={this.handleTaskChange}
              onCreateOption={(inputValue) => this.handleCreateOption('Task', inputValue)}
              options={options.Task} // Set options to the fetched options
            />
            {errors.Task && <p className="error">{errors.Task}</p>}
          </div>

          {/* Downloads */}
          <div>
            <label>Downloads:</label>
            <input
              type="text"
              name="Downloads"
              value={formData.Downloads}
              onChange={this.handleChange}
            />
            {errors.Downloads && <p className="error">{errors.Downloads}</p>}
          </div>

          <button type="submit">Submit</button>
        </form>
      </div>
    );
  }
}

export default MyForm;
