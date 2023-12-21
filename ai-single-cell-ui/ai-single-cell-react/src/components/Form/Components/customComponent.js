import React, { Component , useEffect} from 'react';
import Select from 'react-select';
import CreatableSelect from 'react-select/creatable';
import { SERVER_URL } from '../../../constants/declarations';
import './MyForm.css';
import axios from 'axios';
import { getCookie, isUserAuth } from '../../../utils/utilFunctions';
import { useNavigate } from 'react-router-dom';
import AccessDenied from '../../AccessDeniedPage';

class MyForm extends Component {
  constructor(props) {
    super(props);
    const {
      taskData,
    } = props;

    this.state = {
      formData: taskData.metadata.formData,
      errors: {},
      isLoading: false,
      options: taskData.metadata.options,
      newOptions: taskData.metadata.newOptions,
      message: '',
      hasMessage: false,
      isAdmin: false,
      username: ''
    };  
  }

  componentDidMount() {
    const {taskData} = this.props;

    // Make an API call to get the default options for all fields
    let jwtToken = getCookie('jwtToken');
    if (!jwtToken) {
      // Navigate to the login page using window.location.href
      window.location.href = '/routing';
    } else {
      // If the token exists, verify authenticity
      isUserAuth(jwtToken).then((authData) => {
        if (authData.isAdmin) {
          this.setState({
            isAdmin: authData.isAdmin ,
            username: authData.username, // Set hasMessage to true when a message is set
          });
          if(taskData.metadata.status !== "completed") {
            this.fetchDefaultOptions();
          }
        }
        else {
          console.warn("Unauthorized - you must be an admin to access this page");
          window.location.href = '/accessDenied';
        }
      })
    }
  }

  async fetchDefaultOptions() {
    try {
      const response = await fetch(`${SERVER_URL}/mongoDB/api/options`);
      if (!response.ok) {
        console.error('Error fetching default options');
        return;
      }
      const data = await response.json();
  
      const options = {};
      const fieldNames = [
        'Task', 'Species', 'Sample Type', 'Anatomical Entity',
        'Organ Part', 'Model Organ', 'Selected Cell Types', 'Library Construction Method',
        'Nucleic Acid Source', 'Disease Status (Specimen)', 'Disease Status (Donor)',
        'Development Stage', 'Cell Count Estimate', 'Source',
      ];
  
      fieldNames.forEach(fieldName => {
        if (data[fieldName]) {
          options[fieldName] = data[fieldName].map(option => ({ value: option.abbreviation, label: option.name }));
        }
      });
  
      this.setState({ options });
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

      // Check if the option has already been created to prevent duplicate calls
      if (!this.optionAlreadyCreated(fieldName, inputValue)) {
        this.addNewOptionToMongoDB(fieldName, inputValue);
      }
      this.setState((prevState) => {
        const newOption = { value: inputValue, label: inputValue };
        const updatedOptions = { ...prevState.options };
        updatedOptions[fieldName] = [...(updatedOptions[fieldName] || []), newOption];
    
        const updatedFormData = {
          ...prevState.formData,
          [fieldName]: newOption,
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

  optionAlreadyCreated = (fieldName, inputValue) => {
    return this.state.newOptions.some(
      (option) => option.field === fieldName && option.name === inputValue
    );
  };

  addNewOptionToMongoDB = (fieldName, optionName) => {
    // Make a POST request to your backend to add the new option to MongoDB
    axios
      .post(`${SERVER_URL}/mongoDB/api/addNewOption`, { 'field':fieldName, 'name':optionName, 'username': this.state.username })
      .then((response) => {
        console.log(`New option "${optionName}" added to MongoDB for field "${fieldName}"`);
      })
      .catch((error) => {
        console.error('Error adding new option to MongoDB:', error);
      });
  };

  continueToTaskBuilder = (e) => {
    e.preventDefault();
    const {setFlow} = this.props;
    setFlow('taskBuilder');
  }

  handleSubmit = (e) => {
    e.preventDefault();
    const { setTaskStatus, setTaskData, setActiveTask, taskData} = this.props;

    const errors = this.validateForm(this.state.formData);
    this.setState({ errors });

    if (Object.keys(errors).length === 0) {
      let formData = this.state.formData;

      // add data to the formData
      formData['Cell Count Estimate'] = taskData.quality_control.qc_results[0]?.metadata?.nCells || 0;

      // construct ID 
      // const task_abbv = formData.Task.value;
      const species = formData.Species.value;
      const tissue = formData['Anatomical Entity'].label;
      const cellCount = formData['Cell Count Estimate'];
      const author = formData['Author'];
      const submissionDate = formData['Submission Date'];
      const year = submissionDate ? new Date(submissionDate).getFullYear().toString() : '';

      // Check if cellCount is greater than 1000
      const useCellCount = cellCount && parseInt(cellCount) > 1000;

      const constructedID = `${species}-${tissue}${useCellCount ? `-${cellCount}` : ''}-${author}-${year}`;

      formData.Id = constructedID;

      //Add genes and cells
      formData.Cells = JSON.stringify(taskData.quality_control.qc_results[0]?.metadata?.cells);
      formData.Genes = JSON.stringify(taskData.quality_control.qc_results[0]?.metadata?.genes);

      formData['QC_Plots'] = {
        "scatter_plot": taskData.quality_control.qc_results[0]?.scatter_plot,
        "umap_plot": taskData.quality_control.qc_results[0]?.umap_plot,
        "violin_plot": taskData.quality_control.qc_results[0]?.violin_plot,
        // "highest_expr_genes_plot": taskData.quality_control.qc_results[0]?.highest_expr_genes_plot
      }

      axios.post(`${SERVER_URL}/mongoDB/api/submitDatasetMetadata`, formData)
      .then(response => {
        console.log('Form data submitted successfully');
        this.setState({
          message: 'Dataset created Successfully!',
          hasMessage: true, // Set hasMessage to true when a message is set
        });

        setTaskData((prevTaskData) => ({
          ...prevTaskData,
          metadata: {
            ...prevTaskData.metadata,
            formData: formData,
            taskOptions: this.state.options["Task"],
            options: this.state.options,
            newOptions: this.state.newOptions,
            status: "completed"
          },
        }));
  
          // After Task 1 is successfully completed, update the task status
        setTaskStatus((prevTaskStatus) => ({
          ...prevTaskStatus,
          4: true, // Mark Task 4 as completed
        }));
      })
      .catch(error => {
        console.error('Error submitting form data:', error.response.data.error);
        this.setState({
          message: 'Error submitting form data: ' + error.response.data.error,
          hasMessage: true, // Set hasMessage to true when a message is set
        });
      });

      //The current task is finished, so make the next task active
      // setActiveTask(5);
    } else {
      this.setState({
        message: 'Please fill all the required fields to submit',
        hasMessage: true, // Set hasMessage to true when a message is set
      });
    }
  };

    // Add this method to clear the message and set hasMessage to false
    clearMessageAfterTimeout = () => {
      if (this.state.hasMessage) {
        setTimeout(() => {
          this.setState({
            message: '',
            hasMessage: false,
          });
        }, 5000);
      }
    };

  validateForm(formData) {
    const errors = {};
    if (!formData.Dataset) {
      errors.Dataset = 'Dataset is required';
    }

    if (!formData.Downloads) {
      errors.Downloads = 'Downloads is required';
    }

    if (!formData['Submission Date'] && formData['Submission Date'] === '') {
      errors['Submission Date'] = 'Submission Date is required';
    }

    if (!formData.Author) {
      errors.Author = 'Author is required';
    }

    if (!formData.Species || (formData.Species && formData.Species.value === '')) {
      errors.Species = 'Species is required';
    }

    if (!formData['Organ Part'] || (formData['Organ Part'] && formData['Organ Part'].value === '')) {
      errors['Organ Part'] = 'Organ Part is required';
    }
    if (!formData['Anatomical Entity'] || (formData['Anatomical Entity'] && formData['Anatomical Entity'].value === '')) {
      errors['Anatomical Entity'] = 'Anatomical Entity is required';
    }
    if (!formData['Selected Cell Types'] || (formData['Selected Cell Types'] && formData['Selected Cell Types'].value === '')) {
      errors['Selected Cell Types'] = 'Selected Cell Types is required';
    }
    if (!formData['Disease Status (Specimen)'] || (formData['Disease Status (Specimen)'] && formData['Disease Status (Specimen)'].value === '')) {
      errors['Disease Status (Specimen)'] = 'Disease Status (Specimen) is required';
    }
    if (!formData['Disease Status (Donor)'] || (formData['Disease Status (Donor)'] && formData['Disease Status (Donor)'].value === '')) {
      errors['Disease Status (Donor)'] = 'Disease Status (Donor) is required';
    }
    // if (!formData['Cell Count Estimate'] || (formData['Cell Count Estimate'] && formData['Cell Count Estimate'].value === '')) {
    //   errors['Cell Count Estimate'] = 'Cell Count Estimate is required';
    // }

    const hasError = Object.keys(errors).length > 0;

    return errors;
  }

  render() {
    if (this.state.hasMessage) {
      this.clearMessageAfterTimeout();
    }
    const { formData, errors, isLoading, options, hasMessage, message, isAdmin } = this.state;

    const {setActiveTask, activeTask, taskData } = this.props;
    
    // If isAdmin is false, render nothing
    if (!isAdmin) {
      return null;
    }

    return (
      <div>
        {taskData.metadata.status === 'completed' ? (
          <div>
            <p>Continue to taskBuilder</p>
            <button className="btn btn-info button" onClick={this.continueToTaskBuilder}>Continue</button>
          </div>
        ) : (
      <div className="my-form-container">
        {hasMessage && (
        <div className='message-box' style={{ backgroundColor: '#bdf0c0' }}>
          <div style={{ textAlign: 'center' }}>
            <p>{message}</p>
          </div>
        </div>)}

        <div>
        <h2 className="form-title">My Form</h2>
        <form onSubmit={this.handleSubmit} className="form">
          {/* Dataset */}
          <div className="form-field">
            <label className="form-label">Dataset:</label>
            <input
              type="text"
              name="Dataset"
              required
              value={formData.Dataset}
              onChange={this.handleChange}
              className={`form-input ${errors.Dataset ? 'error' : ''}`}
            />
            {errors.Dataset && <div className="error-tooltip">{errors.Dataset}</div>}
          </div>

          {/* Task (CreatableSelect) */}
          {/* <div className="form-field">
            <label className="form-label">Task:</label>
            <Select
              name="Task"
              value={formData.Task}
              isClearable
              isSearchable
              isLoading={isLoading}
              isCreatable={false} // Set isCreatable to false to prevent creating new options
              onChange={(selectedOption) => this.handleSelectChange('Task', selectedOption)} // Use handleSelectChange            
              options={options.Task} // Set options to the fetched options
              className={`form-input ${errors.Task ? 'error' : ''}`}
            />
            {errors.Task && <div className="error-tooltip">{errors.Task}</div>}
          </div> */}

          {/* Downloads */}
          <div className="form-field">
            <label className="form-label">Downloads:</label>
            <input
              type="text"
              required
              name="Downloads"
              value={formData.Downloads}
              onChange={this.handleChange}
              className={`form-input ${errors.Downloads ? 'error' : ''}`}
            />
            {errors.Downloads && <div className="error-tooltip">{errors.Downloads}</div>}

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
          </div>

          {/* Author (CreatableSelect) */}
          {/* <div className="form-field">
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
              className={`form-input ${errors.Author ? 'error' : ''}`}
            />
            {errors.Author && <div className="error-tooltip">{errors.Author}</div>}
          </div> */}

          <div className="form-field">
            <label className="form-label">Author:</label>
            <input
              type="text"
              name="Author"
              required
              value={formData.Author}
              onChange={this.handleChange}
              className={`form-input ${errors.Author ? 'error' : ''}`}
            />
            {errors.Author && <div className="error-tooltip">{errors.Author}</div>}
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
          </div>

          <div className="form-field">
            <label className="form-label">Abstract:</label>
            <textarea
              name="Abstract"
              value={formData.Abstract}
              onChange={this.handleChange}
              className="form-input"
            />
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
          </div>


          {/* Species (CreatableSelect) */}
          <div className="form-field">
            <label className="form-label">Species:</label>
            <CreatableSelect
              name="Species"
              value={formData.Species}
              isClearable
              isSearchable
              required
              isLoading={isLoading}
              onChange={(selectedOption) => this.handleSelectChange('Species', selectedOption)} // Use handleSelectChange              
              onCreateOption={(inputValue) => this.handleCreateOption('Species', inputValue)}
              options={options.Species} // Set options to the fetched options
              className={`form-input ${errors.Species ? 'error' : ''}`}
            />
            {errors.Species && <div className="error-tooltip">{errors.Species}</div>}
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
              className={`form-input ${errors['Sample Type'] ? 'error' : ''}`}
            />
            {errors['Sample Type'] && <div className="error-tooltip">{errors['Sample Type']}</div>}
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
              className={`form-input ${errors['Anatomical Entity'] ? 'error' : ''}`}
            />
            {errors['Anatomical Entity'] && <div className="error-tooltip">{errors['Anatomical Entity']}</div>}
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
              className={`form-input ${errors['Organ Part'] ? 'error' : ''}`}
            />
            {errors['Organ Part'] && <div className="error-tooltip">{errors['Organ Part']}</div>}
          </div>

          {/* "Model Organ" (CreatableSelect) */}
          <div className="form-field">
            <label className="form-label">Model Organ:</label>
            <CreatableSelect
              name="Model Organ"
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
              className={`form-input ${errors['Selected Cell Types'] ? 'error' : ''}`}
            />
            {errors['Selected Cell Types'] && <div className="error-tooltip">{errors['Selected Cell Types']}</div>}
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
          </div>


          <div className="form-field">
            <label className="form-label">Paired End:</label>
            <div>
              <label>
                <input
                  type="radio"
                  name="Paired End"
                  value= "true"
                  checked={formData["Paired End"] === "true"}
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
                  checked={formData["Paired End"] === "false"}
                  onChange={this.handleChange}
                  className="form-input"
                />
                False
              </label>
            </div>
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
              className={`form-input ${errors['Disease Status (Specimen)'] ? 'error' : ''}`}
            />
            {errors['Disease Status (Specimen)'] && <div className="error-tooltip">{errors['Disease Status (Specimen)']}</div>}
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
              className={`form-input ${errors['Disease Status (Donor)'] ? 'error' : ''}`}
            />
            {errors['Disease Status (Donor)'] && <div className="error-tooltip">{errors['Disease Status (Donor)']}</div>}
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
          </div>


          {/* "Cell Count Estimate" (CreatableSelect) */}
          {/* <div className="form-field">
            <label className="form-label">Cell Count Estimate:</label>
            <input
              type="number"
              name="Cell Count Estimate"
              value={formData['Cell Count Estimate']}
              onChange={this.handleChange}
              className="form-input"
            />
            {errors['Cell Count Estimate'] && <div className="error-tooltip">{errors['Cell Count Estimate']}</div>}
          </div> */}

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
              required
              name="Submission Date"
              value={formData["Submission Date"]}
              onChange={this.handleChange}
              className={`form-input ${errors['Submission Date'] ? 'error' : ''}`}
            />
            {errors['Submission Date'] && <div className="error-tooltip">{errors['Submission Date']}</div>}
          </div>

          <div className='navigation-buttons'>
            <div className="previous">
              <button type="submit" className="btn btn-info button" onClick={() => setActiveTask(activeTask - 1)}>
                Previous
              </button>
            </div>
            <div className="next-upon-success">
              <button type="submit" className="btn btn-info button">
                Submit
              </button>
            </div>
          </div>

        </form>
        </div>
      </div>
        )}
      </div>
    );
  }
}

export default MyForm;
