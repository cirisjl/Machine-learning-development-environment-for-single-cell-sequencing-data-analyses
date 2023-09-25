import React, { Component } from 'react';
import CreatableSelect from 'react-select/creatable';
import { SERVER_URL } from '../../../constants/declarations';
import './MyForm.css';
import axios from 'axios';


class MyForm extends Component {
  constructor(props) {
    super(props);
    this.state = {
      formData: {
        Dataset: '',
        Task: null,
        Downloads: '',
        Title: '',
        Author: null,
        'Reference (paper)':'',
        Abstract: '',
        DOI: '',
        Species: null,
        'Sample Type': null,
        'Anatomical Entity': null,
        'Organ Part': null,
        'Model Organ': null,
        'Selected Cell Types': null,
        'Library Construction Method': null,
        'Nucleic Acid Source': null,
        'Paired End': false,
        'Analysis Protocol': '',
        'Disease Status (Specimen)': null,
        'Disease Status (Donor)': null,
        'Development Stage': null,
        'Donor Count': 0,
        'Cell Count Estimate': 0,
        'Source': null,
        'Source Key': '',
        'Submission Date': 'YYYY-MM-DD', // Set your initial date placeholder here     
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
        'Source': []
      },
      newOptions: [],
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
  
      const updatedFormData = {
        ...prevState.formData,
        [fieldName]: inputValue,
      };
  
      const updatedNewOptions = [
        ...prevState.newOptions,
        { field: fieldName, name: inputValue },
      ];
  
      return {
        options: updatedOptions,
        formData: updatedFormData,
        newOptions: updatedNewOptions,
      };
    });
  };
  

  storeNewOptions = async () => {
    try {
      const { newOptions } = this.state; // Get the new options from your component state
  
      // Assuming you have a variable newOptions containing the new options
      const response = await axios.post(`${SERVER_URL}/mongoDB/api/storeNewOptions`, { newOptions });
  
      if (response.status === 200) {
        console.log('New options stored successfully');
        // If needed, you can update your component state to indicate that new options have been stored.
      } else {
        console.error('Error storing new options');
        // Handle the error as needed
      }
    } catch (error) {
      console.error('Error storing new options:', error);
      // Handle the error as needed
    }
  };
  

  handleSubmit = (e) => {
    e.preventDefault();
    const errors = this.validateForm(this.state.formData);
    this.setState({ errors });

    if (Object.keys(errors).length === 0) {
      const formData = this.state.formData;
      console.log(formData);

      // Assuming you have the form data in a variable named 'formData'
      axios.post(`${SERVER_URL}/mongoDB/api/submitDatasetMetadata`, formData)
        .then(response => {
          console.log('Form data submitted successfully');
          // Insert the new options here by making additional API requests
          this.storeNewOptions();
        })
        .catch(error => {
          console.error('Error submitting form data:', error);
        });
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
      <div className="my-form-container">
        <h2 className="form-title">My Form</h2>
        <form onSubmit={this.handleSubmit} className="form">
          {/* Dataset */}
          <div className="form-field">
            <label className="form-label">Dataset:</label>
            <input
              type="text"
              name="Dataset"
              value={formData.Dataset}
              onChange={this.handleChange}
              className="form-input"
            />
            {errors.Dataset && <p className="error">{errors.Dataset}</p>}
          </div>

          {/* Task (CreatableSelect) */}
          <div className="form-field">
            <label className="form-label">Task:</label>
            <CreatableSelect
              name="Task"
              value={formData.Task}
              isClearable
              isSearchable
              isLoading={isLoading}
              onChange={(selectedOption) => this.handleSelectChange('Task', selectedOption)} // Use handleSelectChange            
              onCreateOption={(inputValue) => this.handleCreateOption('Task', inputValue)}
              options={options.Task} // Set options to the fetched options
              className="form-input"
            />
            {errors.Task && <p className="error">{errors.Task}</p>}
          </div>

          {/* Downloads */}
          <div className="form-field">
            <label className="form-label">Downloads:</label>
            <input
              type="text"
              name="Downloads"
              value={formData.Downloads}
              onChange={this.handleChange}
              className="form-input"
            />
            {errors.Downloads && <p className="error">{errors.Downloads}</p>}
          </div>

          <div className="form-field">
            <label className="form-label">Title:</label>
            <input
              type="text"
              name="Title"
              value={formData.Title}
              onChange={this.handleChange}
              className="form-input"
            />
            {errors.Title && <p className="error">{errors.Title}</p>}
          </div>

          {/* Author (CreatableSelect) */}
          <div className="form-field">
            <label className="form-label">Author:</label>
            <CreatableSelect
              name="Author"
              value={formData.Author}
              isClearable
              isSearchable
              isLoading={isLoading}
              onChange={(selectedOption) => this.handleSelectChange('Author', selectedOption)} // Use handleSelectChange              
              onCreateOption={(inputValue) => this.handleCreateOption('Author', inputValue)}
              options={options.Author} // Set options to the fetched options
              className="form-input"
            />
            {errors.Author && <p className="error">{errors.Author}</p>}
          </div>

          <div className="form-field">
            <label className="form-label">Reference (paper):</label>
            <input
              type="text"
              name="Reference (paper)"
              value={formData['Reference (paper)']}
              onChange={this.handleChange}
              className="form-input"
            />
            {errors['Reference (paper)'] && <p className="error">{errors['Reference (paper)']}</p>}
          </div>

          <div className="form-field">
            <label className="form-label">Abstract:</label>
            <textarea
              name="Abstract"
              value={formData.Abstract}
              onChange={this.handleChange}
              className="form-input"
            />
            {errors.Abstract && <p className="error">{errors.Abstract}</p>}
          </div>

          {/* DOI */}
          <div className="form-field">
            <label className="form-label">DOI:</label>
            <input
              type="text"
              name="DOI"
              value={formData.DOI}
              onChange={this.handleChange}
              placeholder="http://"
              className="form-input"
            />
            {errors.DOI && <p className="error">{errors.DOI}</p>}
          </div>


          {/* Species (CreatableSelect) */}
          <div className="form-field">
            <label className="form-label">Species:</label>
            <CreatableSelect
              name="Species"
              value={formData.Species}
              isClearable
              isSearchable
              isLoading={isLoading}
              onChange={(selectedOption) => this.handleSelectChange('Species', selectedOption)} // Use handleSelectChange              
              onCreateOption={(inputValue) => this.handleCreateOption('Species', inputValue)}
              options={options.Species} // Set options to the fetched options
              className="form-input"
            />
            {errors.Species && <p className="error">{errors.Species}</p>}
          </div>

          {/* "Sample Type" (CreatableSelect) */}
          <div className="form-field">
            <label className="form-label">Sample Type:</label>
            <CreatableSelect
              name="Sample Type"
              value={formData['Sample Type']}
              isClearable
              isSearchable
              isLoading={isLoading}
              onChange={(selectedOption) => this.handleSelectChange('Sample Type', selectedOption)} // Use handleSelectChange             
               onCreateOption={(inputValue) => this.handleCreateOption('Sample Type', inputValue)}
              options={options['Sample Type']} // Set options to the fetched options
              className="form-input"
            />
            {errors['Sample Type'] && <p className="error">{errors['Sample Type']}</p>}
          </div>


          {/* "Anatomical Entity" (CreatableSelect) */}
          <div className="form-field">
            <label className="form-label">Anatomical Entity:</label>
            <CreatableSelect
              name="Anatomical Entity"
              value={formData['Anatomical Entity']}
              isClearable
              isSearchable
              isLoading={isLoading}
              onChange={(selectedOption) => this.handleSelectChange('Anatomical Entity', selectedOption)} // Use handleSelectChange              
              onCreateOption={(inputValue) => this.handleCreateOption('Anatomical Entity', inputValue)}
              options={options['Anatomical Entity']} // Set options to the fetched options
              className="form-input"
            />
            {errors['Anatomical Entity'] && <p className="error">{errors['Anatomical Entity']}</p>}
          </div>

          {/* "Organ Part" (CreatableSelect) */}
          <div className="form-field">
            <label className="form-label">Organ Part:</label>
            <CreatableSelect
              name="Organ Part"
              value={formData['Organ Part']}
              isClearable
              isSearchable
              isLoading={isLoading}
              onChange={(selectedOption) => this.handleSelectChange('Organ Part', selectedOption)} // Use handleSelectChange              
              onCreateOption={(inputValue) => this.handleCreateOption('Organ Part', inputValue)}
              options={options['Organ Part']} // Set options to the fetched options
              className="form-input"
            />
            {errors['Organ Part'] && <p className="error">{errors['Organ Part']}</p>}
          </div>

          {/* "Model Organ" (CreatableSelect) */}
          <div className="form-field">
            <label className="form-label">Model Organ:</label>
            <CreatableSelect
              name="Organ Part"
              value={formData['Model Organ']}
              isClearable
              isSearchable
              isLoading={isLoading}
              onChange={(selectedOption) => this.handleSelectChange('Model Organ', selectedOption)} // Use handleSelectChange              
              onCreateOption={(inputValue) => this.handleCreateOption('Model Organ', inputValue)}
              options={options['Model Organ']} // Set options to the fetched options
              className="form-input"
            />
            {errors['Model Organ'] && <p className="error">{errors['Model Organ']}</p>}
          </div>

          {/* "Selected Cell Types" (CreatableSelect) */}
          <div className="form-field">
            <label className="form-label">Selected Cell Types:</label>
            <CreatableSelect
              name="Selected Cell Types"
              value={formData['Selected Cell Types']}
              isClearable
              isSearchable
              isLoading={isLoading}
              onChange={(selectedOption) => this.handleSelectChange('Selected Cell Types', selectedOption)} // Use handleSelectChange              
              onCreateOption={(inputValue) => this.handleCreateOption('Selected Cell Types', inputValue)}
              options={options['Selected Cell Types']} // Set options to the fetched options
              className="form-input"
            />
            {errors['Selected Cell Types'] && <p className="error">{errors['Selected Cell Types']}</p>}
          </div>



          {/* "Library Construction Method" (CreatableSelect) */}
          <div className="form-field">
            <label className="form-label">Library Construction Method:</label>
            <CreatableSelect
              name="Library Construction Method"
              value={formData['Library Construction Method']}
              isClearable
              isSearchable
              isLoading={isLoading}
              onChange={(selectedOption) => this.handleSelectChange('Library Construction Method', selectedOption)} // Use handleSelectChange              
              onCreateOption={(inputValue) => this.handleCreateOption('Library Construction Method', inputValue)}
              options={options['Library Construction Method']} // Set options to the fetched options
              className="form-input"
            />
            {errors['Library Construction Method'] && <p className="error">{errors['Library Construction Method']}</p>}
          </div>


          {/* "Nucleic Acid Source" (CreatableSelect) */}
          <div className="form-field">
            <label className="form-label">Nucleic Acid Source:</label>
            <CreatableSelect
              name="Nucleic Acid Source"
              value={formData['Nucleic Acid Source']}
              isClearable
              isSearchable
              isLoading={isLoading}
              onChange={(selectedOption) => this.handleSelectChange('Nucleic Acid Source', selectedOption)} // Use handleSelectChange              
              onCreateOption={(inputValue) => this.handleCreateOption('Nucleic Acid Source', inputValue)}
              options={options['Nucleic Acid Source']} // Set options to the fetched options
              className="form-input"
            />
            {errors['Nucleic Acid Source'] && <p className="error">{errors['Nucleic Acid Source']}</p>}
          </div>


          <div className="form-field">
            <label className="form-label">Paired End:</label>
            <div>
              <label>
                <input
                  type="radio"
                  name="Paired End"
                  value="true"
                  checked={formData["Paired End"] === true}
                  onChange={this.handleChange}
                  className="form-input"
                />
                True
              </label>
              <label className="form-label">
                <input
                  type="radio"
                  name="Paired End"
                  value="false"
                  checked={formData["Paired End"] === false}
                  onChange={this.handleChange}
                  className="form-input"
                />
                False
              </label>
            </div>
            {errors["Paired End"] && <p className="error">{errors["Paired End"]}</p>}
          </div>



          
          <div className="form-field">
            <label className="form-label">Analysis Protocol:</label>
            <input
              type="text"
              name="Analysis Protocol"
              value={formData['Analysis Protocol']}
              onChange={this.handleChange}
              className="form-input"
            />
            {errors['Analysis Protocol'] && <p className="error">{errors['Analysis Protocol']}</p>}
          </div>


          
          {/* "Disease Status (Specimen)" (CreatableSelect) */}
          <div className="form-field">
            <label className="form-label">Disease Status (Specimen):</label>
            <CreatableSelect
              name="Disease Status (Specimen)"
              value={formData['Disease Status (Specimen)']}
              isClearable
              isSearchable
              isLoading={isLoading}
              onChange={(selectedOption) => this.handleSelectChange('Disease Status (Specimen)', selectedOption)} // Use handleSelectChange              
              onCreateOption={(inputValue) => this.handleCreateOption('Disease Status (Specimen)', inputValue)}
              options={options['Disease Status (Specimen)']} // Set options to the fetched options
              className="form-input"
            />
            {errors['Disease Status (Specimen)'] && <p className="error">{errors['Disease Status (Specimen)']}</p>}
          </div>


          {/* "Disease Status (Donor)" (CreatableSelect) */}
          <div className="form-field">
            <label className="form-label">Disease Status (Donor):</label>
            <CreatableSelect
              name="Disease Status (Donor)"
              value={formData['Disease Status (Donor)']}
              isClearable
              isSearchable
              isLoading={isLoading}
              onChange={(selectedOption) => this.handleSelectChange('Disease Status (Donor)', selectedOption)} // Use handleSelectChange              
              onCreateOption={(inputValue) => this.handleCreateOption('Disease Status (Donor)', inputValue)}
              options={options['Disease Status (Donor)']} // Set options to the fetched options
              className="form-input"
            />
            {errors['Disease Status (Donor)'] && <p className="error">{errors['Disease Status (Donor)']}</p>}
          </div>


          
          {/* "Development Stage" (CreatableSelect) */}
          <div className="form-field">
            <label className="form-label">Development Stage:</label>
            <CreatableSelect
              name="Development Stage"
              value={formData['Development Stage']}
              isClearable
              isSearchable
              isLoading={isLoading}
              onChange={(selectedOption) => this.handleSelectChange('Development Stage', selectedOption)} // Use handleSelectChange              
              onCreateOption={(inputValue) => this.handleCreateOption('Development Stage', inputValue)}
              options={options['Development Stage']} // Set options to the fetched options
              className="form-input"
            />
            {errors['Development Stage'] && <p className="error">{errors['Development Stage']}</p>}
          </div>

          <div className="form-field">
            <label className="form-label">Donor Count:</label>
            <input
              type="number"
              name="Donor Count"
              value={formData["Donor Count"]}
              onChange={this.handleChange}
              className="form-input"
            />
            {errors["Donor Count"] && <p className="error">{errors["Donor Count"]}</p>}
          </div>


          {/* "Cell Count Estimate" (CreatableSelect) */}
          <div className="form-field">
            <label className="form-label">Cell Count Estimate:</label>
            <CreatableSelect
              name="Cell Count Estimate"
              value={formData['Cell Count Estimate']}
              isClearable
              isSearchable
              isLoading={isLoading}
              onChange={(selectedOption) => this.handleSelectChange('Cell Count Estimate', selectedOption)} // Use handleSelectChange              
              onCreateOption={(inputValue) => this.handleCreateOption('Cell Count Estimate', inputValue)}
              options={options['Cell Count Estimate']} // Set options to the fetched options
              className="form-input"
            />
            {errors['Cell Count Estimate'] && <p className="error">{errors['Cell Count Estimate']}</p>}
          </div>

          {/* "Source" (CreatableSelect) */}
          <div className="form-field">
            <label className="form-label">Source:</label>
            <CreatableSelect
              name="Source"
              value={formData['Source']}
              isClearable
              isSearchable
              isLoading={isLoading}
              onChange={(selectedOption) => this.handleSelectChange('Source', selectedOption)} // Use handleSelectChange              
              onCreateOption={(inputValue) => this.handleCreateOption('Source', inputValue)}
              options={options['Source']} // Set options to the fetched options
              className="form-input"
            />
            {errors['Source'] && <p className="error">{errors['Source']}</p>}
          </div>

          
          {/* Source Key */}
          <div className="form-field">
            <label className="form-label">Source Key:</label>
            <input
              type="text"
              name="Source Key"
              value={formData['Source Key']}
              onChange={this.handleChange}
              placeholder="Enter ..."
              className="form-input"
            />
            {errors['Source Key'] && <p className="error">{errors['Source Key']}</p>}
          </div>

          <div className="form-field">
            <label>Submission Date:</label>
            <input
              type="date"
              name="Submission Date"
              value={formData["Submission Date"]}
              onChange={this.handleChange}
              className="form-input"
            />
            {errors["Submission Date"] && <p className="error">{errors["Submission Date"]}</p>}
          </div>


          <button type="submit">Submit</button>
        </form>
      </div>
    );
  }
}

export default MyForm;
