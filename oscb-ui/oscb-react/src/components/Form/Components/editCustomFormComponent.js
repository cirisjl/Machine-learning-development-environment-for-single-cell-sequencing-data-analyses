import React, { useState, useEffect } from 'react';
import CreatableSelect from 'react-select/creatable';
import { useLocation } from 'react-router-dom';


const EditCustomForm = () => {
  const [formData, setFormData] = useState({});
  const [datasetId, setDatasetId] = useState(null);

  const location = useLocation();

  useEffect(() => {
    const queryParams = new URLSearchParams(location.search);
    const id = queryParams.get('datasetId');
    setDatasetId(id);

    if (id) {
      // Fetch the dataset based on datasetId using axios
      axios.post(`${NODE_API_URL}/item/getDatasetInfoWithPreProcessResults`, { datasetId: id })
        .then(response => {
          setDatasetDetails(response.data);
          setMessage(`Successfully fetched details for the dataset ID - ${id}.`);
          setHasMessage(true);
          setIsError(false);
          setLoading(false);
        })
        .catch(error => {
          console.error(`Error fetching dataset details for dataset ID - ${id}:`, error);
          setMessage(`Error fetching dataset details for dataset ID - ${id}.`);
          setHasMessage(true);
          setIsError(true);
          setLoading(false);
        });
    } else {
      setLoading(false);
    }
  }, [location.search]);

  const [errors, setErrors] = useState({});
  const [hasMessage, setHasMessage] = useState(false);
  const [message, setMessage] = useState('');

  // Handle input changes
  const handleChange = (e) => {
    const { name, value } = e.target;
    setFormData({ ...formData, [name]: value });
  };

  // Handle CreatableSelect changes
  const handleSelectChange = (name, selectedOption) => {
    setFormData({ ...formData, [name]: selectedOption });
  };

  const handleMultiSelectChange = (name, selectedOptions) => {
    setFormData({ ...formData, [name]: selectedOptions || [] });
  };

  const handleCreateOption = (name, inputValue) => {
    setFormData({ ...formData, [name]: inputValue });
  };

  // Handle form submission
  const handleSubmit = (e) => {
    e.preventDefault();
    // Perform validation here
    let formErrors = {};
    if (!formData.Title) formErrors.Title = 'Title is required';
    if (!formData.Author) formErrors.Author = 'Author is required';
    if (!formData.Species) formErrors.Species = 'Species is required';
    if (!formData['Cell Count Estimate']) formErrors['Cell Count Estimate'] = 'Cell Count Estimate is required';
    if (!formData['Organ Part']) formErrors['Organ Part'] = 'Organ Part is required';
    if (!formData['Selected Cell Types'].length) formErrors['Selected Cell Types'] = 'At least one cell type is required';
    if (!formData['Disease Status (Specimen)'].length) formErrors['Disease Status (Specimen)'] = 'At least one disease status is required';
    if (!formData['Submission Date']) formErrors['Submission Date'] = 'Submission date is required';

    if (Object.keys(formErrors).length > 0) {
      setErrors(formErrors);
    } else {
      // Form is valid, submit the data
      setErrors({});
      setHasMessage(true);
      setMessage('Form submitted successfully!');
      // Reset the form if needed
    }
  };

  return (
    <div className="my-form-container">
      {hasMessage && (
        <div className='message-box' style={{ backgroundColor: '#bdf0c0' }}>
          <div style={{ textAlign: 'center' }}>
            <p>{message}</p>
          </div>
        </div>
      )}

      <div>
        <h2 className="form-title">Metadata</h2>
        <form onSubmit={handleSubmit} className="form">
          <div className="form-field">
            <div>
              <label className="form-label">Title:</label>
              <span className="ui-form-title-message warning"> * required, name of the dataset. </span>
            </div>
            <input
              type="text"
              name="Title"
              required
              value={formData.Title}
              onChange={handleChange}
              className={`form-input ${errors.Title ? 'error' : ''}`}
            />
            {errors.Title && <div className="error-tooltip">{errors.Title}</div>}
          </div>

          <div className="form-field">
            <div>
              <label className="form-label">Author:</label>
              <span className="ui-form-title-message warning"> * required </span>
            </div>
            <input
              type="text"
              name="Author"
              required
              value={formData.Author}
              onChange={handleChange}
              className={`form-input ${errors.Author ? 'error' : ''}`}
            />
            {errors.Author && <div className="error-tooltip">{errors.Author}</div>}
          </div>

          {/* Example for a CreatableSelect field */}
          <div className="form-field">
            <div>
              <label className="form-label">Species:</label>
              <span className="ui-form-title-message warning"> * required </span>
            </div>
            <CreatableSelect
              name="Species"
              value={formData.Species}
              isClearable
              isSearchable
              isLoading={isLoading}
              onChange={(selectedOption) => handleSelectChange('Species', selectedOption)}
              onCreateOption={(inputValue) => handleCreateOption('Species', inputValue)}
              options={options.Species}
              className={`form-input ${errors.Species ? 'error' : ''}`}
            />
            {errors.Species && <div className="error-tooltip">{errors.Species}</div>}
          </div>

          <div className='navigation-buttons'>
            <div className="previous">
              <button type="button" className="btn btn-info button" onClick={() => setActiveTask(activeTask - 1)}>
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
  );
};

export default EditCustomForm;
