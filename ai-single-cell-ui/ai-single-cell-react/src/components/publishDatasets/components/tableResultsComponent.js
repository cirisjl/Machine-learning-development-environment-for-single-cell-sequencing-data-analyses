import { faEdit, faEye } from '@fortawesome/free-solid-svg-icons';
import { FontAwesomeIcon } from '@fortawesome/react-fontawesome';
import React, {useState,useEffect} from 'react';
import { useTable, useRowSelect } from 'react-table';
import Button from '@material-ui/core/Button';
import Menu from '@material-ui/core/Menu';
import MenuItem from '@material-ui/core/MenuItem';
import Checkbox from '@material-ui/core/Checkbox';
import FormGroup from '@material-ui/core/FormGroup';
import FormControlLabel from '@material-ui/core/FormControlLabel';
import ListItemText from '@material-ui/core/ListItemText';
import ReactPlotly from './reactPlotly';
import CreatableSelect from 'react-select/creatable';
import {SERVER_URL} from '../../../constants/declarations'
import axios from 'axios';
import '../publishDatasets.css';
import { useNavigate } from 'react-router-dom';





const ResultsTable = ({ data, onSelectDataset, selectedDatasets, multiple, pagination }) => {
  useEffect(() => {
    console.log("data",data);
  });


  const [formData, setFormData] = useState({
    Dataset: '',
    Downloads: '',
    Title: '',
    Author: '',
    'Reference (paper)': '',
    Abstract: '',
    DOI: '',
    Species: '',
    'Sample Type': '',
    'Anatomical Entity': '',
    'Organ Part': '',
    'Model Organ': '',
    'Selected Cell Types': '',
    'Library Construction Method': '',
    'Nucleic Acid Source': '',
    'Paired End': false,
    'Analysis Protocol': '',
    'Disease Status (Specimen)': '',
    'Disease Status (Donor)': '',
    'Development Stage': '',
    'Donor Count': 0,
    'Source': '',
    'Source Key': '',
    'Submission Date': '', // Set your initial date placeholder here
    'Id':''
  });

  const [options,setOptions]=useState({
    Task: [], 
    Author: '',
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
  })

  const [newOption,setNewOption]= useState([]);
  useEffect(()=>{
    console.log(newOption);
  },[newOption])

  const [message, setMessage] = useState('');
const [hasMessage, setHasMessage] = useState(false);
useEffect(() => {
  if (hasMessage) {
    clearMessageAfterTimeout();
  }
}, [hasMessage]); // Run clearMessageAfterTimeout when hasMessage changes

const clearMessageAfterTimeout = () => {
  setTimeout(() => {
    setMessage('');
    setHasMessage(false);
  }, 5000); // 5 seconds timeout
};


  


    const [anchorEl, setAnchorEl] = useState(null);

    // Destructure the pagination object for easier access to its properties
    const { page, pageSize, totalCount } = pagination;

    // Calculate the starting and ending result numbers for the current page
    const startResult = (page - 1) * pageSize + 1;
    const endResult = Math.min(page * pageSize, totalCount); // Ensure not to exceed totalCount
    // const [showView,setShowView]=useState(false);
    const [x,setx]=useState(false);
    const [currentdataset,setCurrentDataset]=useState(null);
    const navigate = useNavigate();


    useEffect(() => {
      console.log("curretn dataset",formData);
    }, [formData]);
  
    const [displayVisualize,setDisplayVisualize] = useState(false);
    const [visualizeDataset, setVisualizeDataset] = useState(null);
    const [activeTab, setActiveTab] = useState(0);

    const handleTabClick = (tabIndex) => {
      setActiveTab(tabIndex);
    };

    
  const handleCloseButtonClick=()=>{
    setDisplayVisualize(false);
  }


    useEffect(() => {
      console.log("visualize dataset",visualizeDataset);
    }, [visualizeDataset]);

   
    useEffect(() => {
      const fetchData = async () => {
        try {
          const response = await fetch(`${SERVER_URL}/mongoDB/api/options`);
          
          if (!response.ok) {
            console.error('Error fetching default options');
            return;
          }
          
          const data = await response.json();
          console.log(data);
    
          const optionss = {};
          const fieldNames = [
            'Task', 'Species', 'Sample Type', 'Anatomical Entity',
            'Organ Part', 'Model Organ', 'Selected Cell Types', 'Library Construction Method',
            'Nucleic Acid Source', 'Disease Status (Specimen)', 'Disease Status (Donor)',
            'Development Stage', 'Cell Count Estimate', 'Source',
          ];
    
          fieldNames.forEach(fieldName => {
            if (data[fieldName]) {
              optionss[fieldName] = data[fieldName].map(option => ({ value: option.abbreviation, label: option.name }));
            }
          });
          console.log(optionss);
          // Update state with fetched options
          setOptions(optionss)

    
        } catch (error) {
          console.error('Error fetching default options:', error);
        }
      };
    
      fetchData(); // Call the fetchData function when the component mounts
    
    }, []); // Empty dependencies array means this effect runs only once after mounting
    
  
    

    const [visibleColumns, setVisibleColumns] = useState({
        'Id': true,
        'TaskId': true,
        'TaskType': true,
        'Title': true,
        'Category': true,
        'Species': true,
        'Organ Part': true,
        'Cell Count Estimate': true,
        'Development Stage': false, // Optional initially not visible
        'Author': false, // Optional initially not visible
        'Submission Date': false, // Optional initially not visible
        'Source': false, // Optional initially not visible
    });

    const handleMenuClick = (event) => {
        setAnchorEl(event.currentTarget);
    };

    const handleMenuClose = () => {
        setAnchorEl(null);
    };

    const toggleColumnVisibility = (column) => {
        console.log("Toggle Column visibility");
        setVisibleColumns(prevVisibleColumns => ({
            ...prevVisibleColumns,
            [column]: !prevVisibleColumns[column],
        }));
        console.log(visibleColumns[column]);
    };

    const resetColumnVisibility = () => {
        setVisibleColumns({
            'Id': true,
            'TaskId': true,
            'TaskType': true,
            'Title': true,
            'Category': true,
            'Species': true,
            'Organ Part': true,
            'Cell Count Estimate': true,
            'Development Stage': false,
            'Author': false,
            'Submission Date': false,
            'Source': false,
        });
    };
    const handleChange = (e, fieldName) => {
      
      const { name, value } = e.target || { name: fieldName, value: e }; // Use fieldName if name is undefined
      console.log(name,value);
      setFormData((prevState) => ({
        ...prevState,
        [name]: value,
      }));
    };

    const handleSelectChange = (fieldName, selectedOption) => {
      setFormData((prevState) => ({
        ...prevState,
        [fieldName]: selectedOption,
      }));
    };
  
    
      
    // const isSelected = datasetId => !!selectedDatasets[datasetId];
    // const isDisabled = () => !multiple && Object.keys(selectedDatasets).length >= 1;

    const handleEdit = (dataset,item) => {
      console.log("dataset----",dataset);
      console.log("item----",item);
      setFormData({
        Dataset:item.Dataset,
        Downloads:item.Downloads,
        Title:item.Title,
        Author:item.Author,
        'Reference (paper)':item['Reference (paper)'],
        Abstract:item.Abstract,
        DOI:item.DOI,
        Species:item.Species,
        'Sample Type':item['Sample Type'],
        'Anatomical Entity': item['Anatomical Entity'],
        'Organ Part': item['Organ Part'],
        'Model Organ': item['Model Organ'],
        'Selected Cell Types': item['Selected Cell Types'],
        'Library Construction Method': item['Library Construction Method'],
        'Nucleic Acid Source': item['Nucleic Acid Source'],
        'Paired End': item['Paired End'],
        'Analysis Protocol': item['Analysis Protocol'],
        'Disease Status (Specimen)': item['Disease Status (Specimen)'],
        'Disease Status (Donor)': item['Disease Status (Donor)'],
        'Development Stage': item['Development Stage'],
        'Donor Count': item['Donor Count'],
        'Source': item['Source'],
        'Source Key': item['Source Key'],
        'Submission Date': item['Submission Date'],
        'Id':dataset

      })
      setx(true);            
    };

    const optionAlreadyCreated = (fieldName, inputValue) => {
      console.log(newOption);
      return newOption.some(
        (option) => option.field === fieldName && option.name === inputValue
      );
    };

    const addNewOptionToMongoDB = (fieldName, optionName) => {
      // Make a POST request to your backend to add the new option to MongoDB
      axios
        .post(`${SERVER_URL}/mongoDB/api/addNewOption`, { 'field':fieldName, 'name':optionName, 'username':'' })
        .then((response) => {
          // console.log(New option "${optionName}" added to MongoDB for field "${fieldName}");
        })
        .catch((error) => {
          console.error('Error adding new option to MongoDB:', error);
        });
    };
     
    const handleCreateOption = (fieldName, inputValue) => {

      // Check if the option has already been created to prevent duplicate calls
      if (! optionAlreadyCreated(fieldName, inputValue)) {
        addNewOptionToMongoDB(fieldName, inputValue);    //Have to implement if required ,the similar api call in CustomComponent.js
      }
      const newOption = { value: inputValue, label: inputValue };
      // Update options state
    const updatedOptions = { ...options };
    updatedOptions[fieldName] = [...(updatedOptions[fieldName] || []), newOption];
    console.log("updated_opti0ons",updatedOptions);
  
  // Update formData state
  const updatedFormData = {
    ...formData,
    [fieldName]: newOption,
  };
  setFormData(updatedFormData);
  console.log(updatedFormData);
  
  setOptions(updatedOptions);
  
  // const updatedNewOptions = [
  //   ...newOption,
  //   { field: fieldName, name: inputValue },
  // ];
  // console.log(updatedNewOptions);
  setNewOption(newOption => [...newOption, { field: fieldName, name: inputValue }]);
  
  // setNewOption(updatedNewOptions);
    
  }; 
   
  
  const handleVisualize = (dataset,id)=>{
    // setVisualizeDataset(dataset);
    // setDisplayVisualize(true);
    // navigate("/handleVisualize/id", { replace: false, state: { newTab: true } });
    const newRoute = `handleVisualize/${id}`
    window.open(`${window.location.origin}/${newRoute}`)



  }

    const handlecloseView=()=>{
      setx(false);
    }
    
  
    const submitHandle = (e) => {
      e.preventDefault();
      console.log("currrent_dataset",currentdataset);
console.log("in handle submit");
console.log(formData);
setVisualizeDataset(formData);

axios.post(`${SERVER_URL}/mongoDB/api/editDatasetMetadata`, formData)
.then(response => {
  console.log('Form data submitted successfully');
  console.log("SERVER_URL",`${SERVER_URL}`)
  setMessage('Dataset Updated Successfully!');
  setHasMessage(true);
})
.catch(error => {
  console.error('Error submitting form data:');
  setMessage('Error submitting form data:');
  setHasMessage(true);

});

    }
    


    const columns = React.useMemo(() => {
        if (data.length === 0) {
            return [];
        }

        const baseColumns = Object.keys(data[0])
        .filter(key => visibleColumns[key])
        .map(key => ({
            Header: key,
            accessor: item => {
                const value = item[key];
                let res = '';
                if (value && typeof value === 'object' && value.label) {
                    res = value.label;
                } else {
                    res = value;
                }
                return(
                    <div 
                        data-title={res} 
                        className="cell-ellipsis"
                        title={res}
                        >
                        {res}
                    </div>
                )
                
            }
        }));

        const actionColumn = {
            id: 'actions',
            Header: 'Actions',
            accessor: item => {
                return(
                    
                <div className="action-buttons">
                    <input
                        type="checkbox"
                        style={{ cursor:'pointer' }}
                        onChange={() => onSelectDataset(item)}
                        checked={!!selectedDatasets[item["Id"]]}
                        // disabled={isDisabled() && !isSelected(item["Id"])} // Disable if multiple is false and a dataset is already selecte

                    />
                    <button
                        onClick={() => handleEdit(item["Id"],item)}
                        className="action-button"
                    >
                        <FontAwesomeIcon icon={faEdit} />
                    </button>
                    <button
                        onClick={() => handleVisualize(item,item["Id"])}
                        className="action-button"
                    >
                        <FontAwesomeIcon icon={faEye} />
                    </button>
                </div>
                );
            }
        };
          
        return [actionColumn, ...baseColumns];
    }, [data, selectedDatasets, visibleColumns]);

    const {
        getTableProps,
        getTableBodyProps,
        headerGroups,
        rows,
        prepareRow,
    } = useTable({ columns, data },useRowSelect);

    return (
    <div>
        <div>
            {/* Dropdown for editing columns */}
            <div className="dropdown">
                <div className='total-results-count'>
                <p>Results {startResult} - {endResult} of {totalCount}</p>
                </div>
                <Button variant="contained" aria-controls="edit-columns-menu" aria-haspopup="true" onClick={handleMenuClick} >
                    Edit Columns
                </Button>
                <Menu
                    id="edit-columns-menu"
                    anchorEl={anchorEl}
                    keepMounted
                    open={Boolean(anchorEl)}
                    onClose={handleMenuClose}
                    style={{height: "400px" }} // Adjusts the vertical position
                    anchorOrigin={{
                        vertical: 'bottom', 
                        horizontal: 'right', 
                      }}
                      transformOrigin={{
                        vertical: 'top', 
                        horizontal: 'right',
                      }}
                      getContentAnchorEl={null} // This will make anchorOrigin work as expected
                >
                    <FormGroup>
                        {Object.keys(visibleColumns).map((column) => (
                           <MenuItem key={column} onClick={(event) => event.stopPropagation()}>
                           <FormControlLabel
                               control={
                                   <Checkbox
                                       checked={visibleColumns[column]}
                                       onChange={() => toggleColumnVisibility(column)}
                                       onClick={(event) => event.stopPropagation()} // Prevent triggering the menu item's onClick
                                   />
                               }
                               label={column}
                               // Remove the onClick here to avoid overriding Checkbox's behavior
                           />
                       </MenuItem>
                        ))}
                        {/* Reset Menu Item */}
                        <MenuItem onClick={resetColumnVisibility}>
                            <ListItemText primary="Reset" />
                        </MenuItem>
                    </FormGroup>
                </Menu>
            </div>

            <table {...getTableProps()} className="table-container">
            <   thead>
                    {headerGroups.map((headerGroup, index) => (
                        <tr {...headerGroup.getHeaderGroupProps()} key={index}>
                            {headerGroup.headers.map((column, colIndex) => (
                                <th {...column.getHeaderProps()} key={colIndex}>{column.render('Header')}</th>
                            ))}
                        </tr>
                    ))}
                </thead>
                <tbody {...getTableBodyProps()}>
                    {rows.map((row, rowIndex) => {
                        prepareRow(row);
                        return (
                            <tr {...row.getRowProps()} key={rowIndex}>
                                {row.cells.map((cell, cellIndex) => {
                                    return <td {...cell.getCellProps()} key={cellIndex}>{cell.value}</td>;
                                })}
                            </tr>
                        );
                    })}
                </tbody>
            </table>
            
        </div>
        <div>

        {x && (
          <div className="dialog-container">
          <div className="dialog1">

          <div className="my-form-container1"style={{ height: '600px'}}>
          {hasMessage && (
        <div className='message-box' style={{ backgroundColor: '#bdf0c0' }}>
          <div style={{ textAlign: 'center' }}>
            <p>{message}</p>
          </div>
        </div>)}
        
        <div>
        <h2 className="form-title">My Form</h2>
        <form  className="form1" >
          {/* Dataset */}
          <div className="form-field1">
            <label className="form-label1">Dataset:</label>
            <input
              type="text"
              name="Dataset"
              value={formData.Dataset}
              onChange={(e) => handleChange(e, 'Dataset')}
              

            //   required
            />
          </div>

          

          {/* Downloads */}
          <div className="form-field1">
            <label className="form-label1">Downloads:</label>
            <input
              type="text"
              required
              name="Downloads"
              value={formData.Downloads}
              onChange={handleChange}

            />

          </div>

          <div className="form-field1">
            <label className="form-label1">Title:</label>
            <input
              type="text"
              name="Title"
              value={formData.Title}
              onChange={handleChange}

              //className="form-input"
            />
          </div>

          <div className="form-field1">
            <label className="form-label1">Author:</label>
            <input
              type="text"
              name="Author"
              value={formData.Author}
              onChange={handleChange}

              required
              
            />
          </div>

          <div className="form-field1">
            <label className="form-label1">Reference (paper):</label>
            <input
              type="text"
              name="Reference (paper)"
              value={formData['Reference (paper)']}
              onChange={handleChange}

              className="form-input"
            />
          </div>

          <div className="form-field1">
            <label className="form-label1">Abstract:</label>
            <textarea
              name="Abstract"
              value={formData.Abstract}
              onChange={handleChange}

              className="form-input"
            />
          </div>

          <div className="form-field1">
            <label className="form-label1">DOI:</label>
            <input
              type="text"
              name="DOI"
              value={formData.DOI}
              onChange={handleChange}

              placeholder="http://"
              className="form-input"
            //   value={currentdataset.}
            />
          </div>


          {/* Species (CreatableSelect) */}
          <div className="form-field1">
            <label className="form-label1">Species:</label>
            <CreatableSelect
              name="Species"
              value={formData.Species}
              isClearable
              isSearchable
              required
              onChange={(selectedOption) => handleSelectChange('Species', selectedOption)} // Use handleSelectChange              
              onCreateOption={(inputValue) => handleCreateOption('Species', inputValue)}
              options={options.Species} // Set options to the fetched options
            />
          </div>

          {/* "Sample Type" (CreatableSelect) */}
          <div className="form-field1">
            <label className="form-label1">Sample Type:</label>
            <CreatableSelect
              name="Sample Type"
              value={formData['Sample Type']}

              isClearable
              isSearchable
            //   isLoading={isLoading}
              onChange={(selectedOption) => handleSelectChange('Sample Type', selectedOption)} // Use handleSelectChange             
               onCreateOption={(inputValue) => handleCreateOption('Sample Type', inputValue)}
              options={options['Sample Type']} // Set options to the fetched options
            />
          </div>


          {/* "Anatomical Entity" (CreatableSelect) */}
          <div className="form-field1">
            <label className="form-label1">Anatomical Entity:</label>
            <CreatableSelect
              name="Anatomical Entity"
              value={formData['Anatomical Entity']}

              isClearable
              isSearchable
            //   isLoading={isLoading}
              onChange={(selectedOption) => handleSelectChange('Anatomical Entity', selectedOption)} // Use handleSelectChange              
              onCreateOption={(inputValue) => handleCreateOption('Anatomical Entity', inputValue)}
              options={options['Anatomical Entity']} // Set options to the fetched options
            />
          </div>

          {/* "Organ Part" (CreatableSelect) */}
          <div className="form-field1">
            <label className="form-label1">Organ Part:</label>
            <CreatableSelect
              name="Organ Part"
              value={formData['Organ Part']}
              isClearable
              isSearchable
            //   isLoading={isLoading}
              onChange={(selectedOption) => handleSelectChange('Organ Part', selectedOption)} // Use handleSelectChange              
              onCreateOption={(inputValue) => handleCreateOption('Organ Part', inputValue)}
              options={options['Organ Part']} // Set options to the fetched options
            />
                 </div>

          {/* "Model Organ" (CreatableSelect) */}
          <div className="form-field1">
            <label className="form-label1">Model Organ:</label>
            <CreatableSelect
              name="Model Organ"
              value={formData['Model Organ']}

              isClearable
              isSearchable
            //   isLoading={isLoading}
              onChange={(selectedOption) => handleSelectChange('Model Organ', selectedOption)} // Use handleSelectChange              
              onCreateOption={(inputValue) => handleCreateOption('Model Organ', inputValue)}
              options={options['Model Organ']} // Set options to the fetched options
              className="form-input"
            />
          </div>

          {/* "Selected Cell Types" (CreatableSelect) */}
          <div className="form-field1">
            <label className="form-label1">Selected Cell Types:</label>
            <CreatableSelect
              name="Selected Cell Types"
              value={formData['Selected Cell Types']}

              isClearable
              isSearchable
            //   isLoading={isLoading}
              onChange={(selectedOption) => handleSelectChange('Selected Cell Types', selectedOption)} // Use handleSelectChange              
              onCreateOption={(inputValue) => handleCreateOption('Selected Cell Types', inputValue)}
              options={options['Selected Cell Types']} // Set options to the fetched options
            />
          </div>



          {/* "Library Construction Method" (CreatableSelect) */}
          <div className="form-field1">
            <label className="form-label1">Library Construction Method:</label>
            <CreatableSelect
              name="Library Construction Method"
              value={formData['Library Construction Method']}

              isClearable
              isSearchable
            //   isLoading={isLoading}
              onChange={(selectedOption) => handleSelectChange('Library Construction Method', selectedOption)} // Use handleSelectChange              
              onCreateOption={(inputValue) => handleCreateOption('Library Construction Method', inputValue)}
              options={options['Library Construction Method']} // Set options to the fetched options
              className="form-input"
            />
          </div>


          {/* "Nucleic Acid Source" (CreatableSelect) */}
          <div className="form-field1">
            <label className="form-label1">Nucleic Acid Source:</label>
            <CreatableSelect
              name="Nucleic Acid Source"
              value={formData['Nucleic Acid Source']}

              isClearable
              isSearchable
            //   isLoading={isLoading}
              onChange={(selectedOption) => handleSelectChange('Nucleic Acid Source', selectedOption)} // Use handleSelectChange              
              onCreateOption={(inputValue) => handleCreateOption('Nucleic Acid Source', inputValue)}
              options={options['Nucleic Acid Source']} // Set options to the fetched options
              className="form-input"
            />
          </div>


          <div className="form-field1">
            <label className="form-label1">Paired End:</label>
            <div>
              <label>
                <input
                  type="radio"
                  name="Paired End"
                  value= "true"
                  checked={formData["Paired End"] === "true"}
                onChange={handleChange}

                  className="form-input"
                />
                True
              </label>
              <label className="form-label1">
                <input
                  type="radio"
                  name="Paired End"
                  value="false"
                  checked={formData["Paired End"] === "false"}
                  className="form-input"
                  onChange={handleChange}

                />
                False
              </label>
            </div>
          </div>

          <div className="form-field1">
            <label className="form-label1">Analysis Protocol:</label>
            <input
              type="text"
              name="Analysis Protocol"
              value={formData['Analysis Protocol']}
            //   onChange={this.handleChange}
            onChange={handleChange}

              className="form-input"
            />
          </div>

          {/* "Disease Status (Specimen)" (CreatableSelect) */}
          <div className="form-field1">
            <label className="form-label1">Disease Status (Specimen):</label>
            <CreatableSelect
              name="Disease Status (Specimen)"
              value={formData['Disease Status (Specimen)']}

              isClearable
              isSearchable
            //   isLoading={isLoading}
              onChange={(selectedOption) => handleSelectChange('Disease Status (Specimen)', selectedOption)} // Use handleSelectChange              
              onCreateOption={(inputValue) => handleCreateOption('Disease Status (Specimen)', inputValue)}
              options={options['Disease Status (Specimen)']} // Set options to the fetched options
            />
            {/* {errors['Disease Status (Specimen)'] && <div className="error-tooltip">{errors['Disease Status (Specimen)']}</div>} */}
          </div>


          {/* "Disease Status (Donor)" (CreatableSelect) */}
          <div className="form-field1">
            <label className="form-label1">Disease Status (Donor):</label>
            <CreatableSelect
              name="Disease Status (Donor)"
              value={formData['Disease Status (Donor)']}

              isClearable
              isSearchable
            //   isLoading={isLoading}
              onChange={(selectedOption) => handleSelectChange('Disease Status (Donor)', selectedOption)} // Use handleSelectChange              
              onCreateOption={(inputValue) => handleCreateOption('Disease Status (Donor)', inputValue)}
              options={options['Disease Status (Donor)']} // Set options to the fetched options
            />
            {/* {errors['Disease Status (Donor)'] && <div className="error-tooltip">{errors['Disease Status (Donor)']}</div>} */}
          </div>

          {/* "Development Stage" (CreatableSelect) */}
          <div className="form-field1">
            <label className="form-label1">Development Stage:</label>
            <CreatableSelect
              name="Development Stage"
              value={formData['Development Stage']}

              isClearable
              isSearchable
            //   isLoading={isLoading}
              onChange={(selectedOption) => handleSelectChange('Development Stage', selectedOption)} // Use handleSelectChange              
              onCreateOption={(inputValue) => handleCreateOption('Development Stage', inputValue)}
              options={options['Development Stage']} // Set options to the fetched options
              className="form-input"
            />
          </div>

          <div className="form-field1">
            <label className="form-label1">Donor Count:</label>
            <input
              type="number"
              name="Donor Count"
              value={formData["Donor Count"]}
            //   onChange={this.handleChange}
            onChange={handleChange}

              className="form-input"
            />
          </div>


         

          {/* "Source" (CreatableSelect) */}
          <div className="form-field1">
            <label className="form-label1">Source:</label>
            <CreatableSelect
              name="Source"
              value={formData['Source']}

              isClearable
              isSearchable
            //   isLoading={isLoading}
              onChange={(selectedOption) => handleSelectChange('Source', selectedOption)} // Use handleSelectChange              
              onCreateOption={(inputValue) => handleCreateOption('Source', inputValue)}
              options={options['Source']} // Set options to the fetched options
              className="form-input"
            />
          </div>

          
          {/* Source Key */}
          <div className="form-field1">
            <label className="form-label1">Source Key:</label>
            <input
              type="text"
              name="Source Key"
              value={formData['Source Key']}
            //   onChange={this.handleChange}
            onChange={handleChange}

              placeholder="Enter ..."
              className="form-input"
            />
            {/* {errors['Source Key'] && <p className="error">{errors['Source Key']}</p>} */}
          </div>

          <div className="form-field1">
            <label>Submission Date:</label>
            <input
              type="date"
              required
              name="Submission Date"
              value={formData["Submission Date"]}
              onChange={handleChange}

            //   onChange={this.handleChange}
            />
            {/* {errors['Submission Date'] && <div className="error-tooltip">{errors['Submission Date']}</div>} */}
          </div>

          <div style={{ display: 'flex', justifyContent: 'center', marginTop: '20px' }}>
  <div style={{ marginRight: '10px' }}>
    <button type="submit" onClick={submitHandle}>Save</button>
  </div>
  <div>
    <button onClick={handlecloseView}>Close</button>
  </div>
</div>





         

         

          

        </form>
        </div>
      </div>
           

           
            </div>
          </div>
        )}
        {displayVisualize && (
          <div className="dialog-container">
          <div className="dialog1">
          <div className="tab">
        <button className={activeTab === 0 ? 'active' : ''} onClick={() => handleTabClick(0)}>Summary</button>
        <button className={activeTab === 1 ? 'active' : ''} onClick={() => handleTabClick(1)}>Explore</button>
        <button className={activeTab === 2 ? 'active' : ''} onClick={() => handleTabClick(2)}>Download</button>
         </div>
      <div className="tab-content">
        {activeTab === 0 && 
        <div className="tab-panel">
           <h2>Title</h2>
            <p>{visualizeDataset.Title}</p>
          <hr />
          <h2>Abstract</h2>
          <p>{visualizeDataset.Abstract}</p>
          <hr />
        <h2>Author</h2>
        <p>{visualizeDataset.Author}</p>
        <hr />
          <h2>Reference (paper)</h2>
          <p>{visualizeDataset['Reference (paper)']}</p>
          <hr />
          <h2>Species</h2>
          <p>{visualizeDataset['Species'].value}</p>
          <hr/>
          <h2>Sample Type</h2>
          <p>{visualizeDataset['Sample Type'].value}</p>
          <hr/>
          <p>{visualizeDataset['Anatomical Entity'].label}</p>
          <hr/>
          <h2>Organ Part</h2>
          <p>{visualizeDataset['Organ Part'].label}</p>
          <hr/>
          <h2>Model Organ</h2>
          <p>{visualizeDataset['Model Organ'].label}</p>
          <br/>
          <h2>Selected Cell Types</h2>
          <p>{visualizeDataset['Selected Cell Types'].label}</p>
          <hr/>
          <h2>Library Construction Method</h2>
          <p>{visualizeDataset['Library Construction Method']}</p>
          <hr/>
          <h2>Nucleic Acid Source</h2>
          <p>{visualizeDataset['Nucleic Acid Source']}</p>
          <hr/>
          <h2>Disease Status (Specimen)</h2>
          <p>{visualizeDataset['Disease Status (Specimen)'].label}</p>
          <hr/>
          <h2>Disease Status (Donor)</h2>
          <p>{visualizeDataset['Disease Status (Donor)'].label}</p>
          <hr/>
          <h2>Development Stage</h2>
          <p>{visualizeDataset['Development Stage']}</p>
          <hr/>
          <h2>Donor Count</h2>
          <p>{visualizeDataset['Donor Count']}</p>
          <hr/>
          <h2>Source</h2>
          <p>{visualizeDataset['Source']}</p>
          
          <button className="close-button" onClick={() => handleCloseButtonClick()}>Close</button>

          </div>}
        {activeTab === 1 && <div className="tab-panel">
          <h2>Scatter Plot</h2>
          <ReactPlotly plot_data={visualizeDataset.QC_Plots.scatter_plot} />

          <h2>UMAP Plot</h2>
          <ReactPlotly plot_data={visualizeDataset.QC_Plots.umap_plot} />


          <h2>Violin Plot</h2>
          <ReactPlotly plot_data={visualizeDataset.QC_Plots.violin_plot} />


          <button className="close-button" onClick={() => handleCloseButtonClick()}>Close</button>


          </div>}
        {activeTab === 2 && <div className="tab-panel">Content for Tab 3
          <button className="close-button" onClick={() => handleCloseButtonClick()}>Close</button>
</div>}
      </div>

            </div>
            </div>

        )}
        </div>
        </div>
    );
};

export default ResultsTable;
