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
        Task: [], // Initialize Task options as an empty array
      },
    };
  }

  componentDidMount() {
    // Make an API call to get the default options for all fields
    this.fetchDefaultOptions();
  }

  async fetchDefaultOptions() {
    try {
      const response = await fetch(`${SERVER_URL}/mongoDB/api/options`);
      if (response.ok) {
        const data = await response.json();

        this.setState({
          options: {
            Task: data.Task.map((option) => ({ value: option, label: option })),
          },
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
