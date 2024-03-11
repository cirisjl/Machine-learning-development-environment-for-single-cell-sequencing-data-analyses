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
import ReactDOM from 'react-dom';
import { createRoot } from 'react-dom/client';
import CreatableSelect from 'react-select/creatable';
// import { SERVER_URL } from '../../../constants/declarations';
import {SERVER_URL} from '../../../constants/declarations'
import axios from 'axios';




const ResultsTable = ({ data, onSelectDataset, selectedDatasets, multiple, pagination }) => {
  useEffect(() => {
    console.log("data",data);
  }, [data]);


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
  });




  


    const [anchorEl, setAnchorEl] = useState(null);

    // Destructure the pagination object for easier access to its properties
    const { page, pageSize, totalCount } = pagination;

    // Calculate the starting and ending result numbers for the current page
    const startResult = (page - 1) * pageSize + 1;
    const endResult = Math.min(page * pageSize, totalCount); // Ensure not to exceed totalCount
    // const [showView,setShowView]=useState(false);
    const [x,setx]=useState(false);
    const [currentdataset,setCurrentDataset]=useState(null);
    useEffect(() => {
      console.log("curretn dataset",formData);
    }, [formData]);
  
    

    const [visibleColumns, setVisibleColumns] = useState({
        'Id': true,
        'TaskId': true,
        'TaskType': true,
        'Title': true,
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
    
      
    const isSelected = datasetId => !!selectedDatasets[datasetId];
    const isDisabled = () => !multiple && Object.keys(selectedDatasets).length >= 1;

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
        'Submission Date': item['Submission Date']

      })
      setx(true);            
    };
    const optionAlreadyCreated = (fieldName, inputValue) => {
      return this.state.newOptions.some(
        (option) => option.field === fieldName && option.name === inputValue
      );
    };
   
  const handleCreateOption = (fieldName, inputValue) => {

    // Check if the option has already been created to prevent duplicate calls
    if (! optionAlreadyCreated(fieldName, inputValue)) {
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
    const handleVisualize = (dataset) => {
        console.log("Visualizing dataset: ", dataset);
        console.log("resulys", data);
        console.log("hghg",data.QC_Plots)
        console.log("check2",data[0].QC_Plots);
        console.log("check2",data[0].QC_Plots.scatter_plot);
        
        // Open a new tab/window
        const newWindow = window.open('', '_blank');
    
        if (newWindow) {
            // Write the necessary HTML structure to the new tab with a spinner
            newWindow.document.write('<html><head><title>Loading...</title>');
            newWindow.document.write('<style>');
            newWindow.document.write(`
                body { display: flex; flex-direction: column; margin: 0; }
                .spinner { border: 5px solid #f3f3f3; border-top: 5px solid #3498db; border-radius: 50%; width: 50px; height: 50px; animation: spin 1s linear infinite; }
                .plot-container { flex: 1; display: flex; flex-direction: column; justify-content: center; align-items: center; height: 100vh; margin-top: 20px; }
                .plot-container h2 { margin-bottom: 20px; }
            `);
            newWindow.document.write('</style></head>');
            newWindow.document.write('<body><div class="spinner"></div></body></html>');
    
            // Debug statement to check if the spinner HTML is written
            console.log("Spinner HTML written to the new tab.");
    
            setTimeout(() => {
                if (data && data.length > 0 && data[0].QC_Plots) {
                    // Debug statement to check the availability of QC_Plots in the data
                    console.log("QC_Plots found in data.");
    
                    // Write the necessary HTML structure to the new tab
                    newWindow.document.write('<html><head><title>Visualization</title></head><body></body></html>');
    
                    // Debug statement to verify the HTML structure is written correctly
                    console.log("HTML structure written to the new tab.");
    
                    // Render your plots in the new window/tab, for example using React components:
                    if (data[0].QC_Plots.umap_plot) {
                        const umapPlotContainer = newWindow.document.createElement('div');
                        umapPlotContainer.className = 'plot-container';
                        umapPlotContainer.innerHTML = '<h2>UMAP Plot</h2>';
                        newWindow.document.body.appendChild(umapPlotContainer);
                        const root = newWindow.document.createElement('div');
                        umapPlotContainer.appendChild(root);
                        createRoot(root).render(<ReactPlotly plot_data={data[0].QC_Plots.umap_plot} />);
                    }
    
                    if (data[0].QC_Plots.scatter_plot) {
                        const scatterPlotContainer = newWindow.document.createElement('div');
                        scatterPlotContainer.className = 'plot-container';
                        scatterPlotContainer.innerHTML = '<h2>Scatter Plot</h2>';
                        newWindow.document.body.appendChild(scatterPlotContainer);
                        const root = newWindow.document.createElement('div');
                        scatterPlotContainer.appendChild(root);
                        createRoot(root).render(<ReactPlotly plot_data={data[0].QC_Plots.scatter_plot} />);
                    }
    
                    if (data[0].QC_Plots.violin_plot) {
                        const violinPlotContainer = newWindow.document.createElement('div');
                        violinPlotContainer.className = 'plot-container';
                        violinPlotContainer.innerHTML = '<h2>Violin Plot</h2>';
                        newWindow.document.body.appendChild(violinPlotContainer);
                        const root = newWindow.document.createElement('div');
                        violinPlotContainer.appendChild(root);
                        createRoot(root).render(<ReactPlotly plot_data={data[0].QC_Plots.violin_plot} />);
                    }
                    const spinner = newWindow.document.querySelector(".spinner");
                spinner.parentNode.removeChild(spinner);
                } else {
                    // Error handling if data or QC_Plots is not found
                    console.error("Data or QC_Plots not found or not in the expected format.");
                }
    
                // For demonstration, let's write a simple message after a delay
                //newWindow.document.body.innerHTML = '<h1>Data Loaded!</h1>';
            }, 2000); // Simulating a delay of 2 seconds (2000 milliseconds)
        } else {
            // Error handling if newWindow is null (window.open failed)
            console.error("Failed to open new window. Ensure popup blockers are disabled.");
        }
    };
    

    const handlecloseView=()=>{
      setx(false);
    }
    
    const submitHandle = (e) => {
      e.preventDefault();
      console.log("currrent_dataset",currentdataset);
      // setFormData( {
      //   Dataset: currentdataset.Dataset,
      //   Downloads: currentdataset.Downloads,
      //   Title: currentdataset.Title,
      //   Author: currentdataset.Author,
      //   'Reference (paper)': currentdataset['Reference (paper)'],
      //   Abstract: currentdataset.Abstract,
      //   DOI: currentdataset.DOI,
      //   Species: currentdataset.Species,
      //   'Sample Type': currentdataset['Sample Type'],
      //   'Anatomical Entity': currentdataset['Anatomical Entity'],
      //   'Organ Part': currentdataset['Organ Part'],
      //   'Model Organ': currentdataset['Model Organ'],
      //   'Selected Cell Types': currentdataset['Selected Cell Types'],
      //   'Library Construction Method': currentdataset['Library Construction Method'],
      //   'Nucleic Acid Source': currentdataset['Nucleic Acid Source'],
      //   'Paired End': currentdataset['Paired End'],
      //   'Analysis Protocol': currentdataset['Analysis Protocol'],
      //   'Disease Status (Specimen)': currentdataset['Disease Status (Specimen)'],
      //   'Disease Status (Donor)': currentdataset['Disease Status (Donor)'],
      //   'Development Stage': currentdataset['Development Stage'],
      //   'Donor Count': currentdataset['Donor Count'],
      //   'Source': currentdataset.Source,
      //   'Source Key': currentdataset['Source Key'],
      //   'Submission Date': currentdataset['Submission Date'],
      // });
      
console.log("in handle submit");
console.log(formData);


axios.post(`${SERVER_URL}/mongoDB/api/editDatasetMetadata`, formData)
.then(response => {
  console.log('Form data submitted successfully');
  this.setState({
    message: 'Dataset created Successfully!',
    hasMessage: true, // Set hasMessage to true when a message is set
  });


 
})
.catch(error => {
  console.error('Error submitting form data:', error.response.data.error);
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
                        disabled={isDisabled() && !isSelected(item["Id"])} // Disable if multiple is false and a dataset is already selecte

                    />
                    <button
                        onClick={() => handleEdit(item["Id"],item)}
                        className="action-button"
                    >
                      
                        <FontAwesomeIcon icon={faEdit} />
                    </button>
                    <button
                        onClick={() => handleVisualize(item["Id"])}
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
          <div  className="dialog-container">
          <div className="dialog">

          <div className="my-form-container">
        
        <div>
        <h2 className="form-title">My Form</h2>
        <form  className="form" >
          {/* Dataset */}
          <div className="form-field">
            <label className="form-label">Dataset:</label>
            <input
              type="text"
              name="Dataset"
              value={formData.Dataset}
              onChange={(e) => handleChange(e, 'Dataset')}
              

            //   required
            />
          </div>

          

          {/* Downloads */}
          <div className="form-field">
            <label className="form-label">Downloads:</label>
            <input
              type="text"
              required
              name="Downloads"
              value={formData.Downloads}
              onChange={handleChange}

            />

          </div>

          <div className="form-field">
            <label className="form-label">Title:</label>
            <input
              type="text"
              name="Title"
              value={formData.Title}
              onChange={handleChange}

              //className="form-input"
            />
          </div>

          <div className="form-field">
            <label className="form-label">Author:</label>
            <input
              type="text"
              name="Author"
              value={formData.Author}
              onChange={handleChange}

              required
              
            />
          </div>

          <div className="form-field">
            <label className="form-label">Reference (paper):</label>
            <input
              type="text"
              name="Reference (paper)"
              value={formData['Reference (paper)']}
              onChange={handleChange}

              className="form-input"
            />
          </div>

          <div className="form-field">
            <label className="form-label">Abstract:</label>
            <textarea
              name="Abstract"
              value={formData.Abstract}
              onChange={handleChange}

              className="form-input"
            />
          </div>

          <div className="form-field">
            <label className="form-label">DOI:</label>
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
          <div className="form-field">
            <label className="form-label">Species:</label>
            <CreatableSelect
              name="Species"
              value={formData.Species}
              onChange={handleChange}

              isClearable
              isSearchable
              required
            //   isLoading={isLoading}
            //   onChange={(selectedOption) => this.handleSelectChange('Species', selectedOption)} // Use handleSelectChange              
              onCreateOption={(inputValue) => handleCreateOption('Species', inputValue)}
            //   options={options.Species} // Set options to the fetched options
            />
          </div>

          {/* "Sample Type" (CreatableSelect) */}
          <div className="form-field">
            <label className="form-label">Sample Type:</label>
            <CreatableSelect
              name="Sample Type"
              value={FormData['Sample Type']}
              onChange={handleChange}

              isClearable
              isSearchable
            //   isLoading={isLoading}
            //   onChange={(selectedOption) => this.handleSelectChange('Sample Type', selectedOption)} // Use handleSelectChange             
               onCreateOption={(inputValue) => handleCreateOption('Sample Type', inputValue)}
            //   options={options['Sample Type']} // Set options to the fetched options
            />
          </div>


          {/* "Anatomical Entity" (CreatableSelect) */}
          <div className="form-field">
            <label className="form-label">Anatomical Entity:</label>
            <CreatableSelect
              name="Anatomical Entity"
              value={formData['Anatomical Entity']}
              onChange={handleChange}

              isClearable
              isSearchable
            //   isLoading={isLoading}
            //   onChange={(selectedOption) => this.handleSelectChange('Anatomical Entity', selectedOption)} // Use handleSelectChange              
            //   onCreateOption={(inputValue) => this.handleCreateOption('Anatomical Entity', inputValue)}
            //   options={options['Anatomical Entity']} // Set options to the fetched options
            />
          </div>

          {/* "Organ Part" (CreatableSelect) */}
          <div className="form-field">
            <label className="form-label">Organ Part:</label>
            <CreatableSelect
              name="Organ Part"
              value={formData['Organ Part']}
              onChange={handleChange}

              isClearable
              isSearchable
            //   isLoading={isLoading}
            //   onChange={(selectedOption) => this.handleSelectChange('Organ Part', selectedOption)} // Use handleSelectChange              
            //   onCreateOption={(inputValue) => this.handleCreateOption('Organ Part', inputValue)}
            //   options={options['Organ Part']} // Set options to the fetched options
            />
                 </div>

          {/* "Model Organ" (CreatableSelect) */}
          <div className="form-field">
            <label className="form-label">Model Organ:</label>
            <CreatableSelect
              name="Model Organ"
              value={formData['Model Organ']}
              onChange={handleChange}

              isClearable
              isSearchable
            //   isLoading={isLoading}
            //   onChange={(selectedOption) => this.handleSelectChange('Model Organ', selectedOption)} // Use handleSelectChange              
            //   onCreateOption={(inputValue) => this.handleCreateOption('Model Organ', inputValue)}
            //   options={options['Model Organ']} // Set options to the fetched options
              className="form-input"
            />
          </div>

          {/* "Selected Cell Types" (CreatableSelect) */}
          <div className="form-field">
            <label className="form-label">Selected Cell Types:</label>
            <CreatableSelect
              name="Selected Cell Types"
              value={FormData['Selected Cell Types']}
              onChange={handleChange}

              isClearable
              isSearchable
            //   isLoading={isLoading}
            //   onChange={(selectedOption) => this.handleSelectChange('Selected Cell Types', selectedOption)} // Use handleSelectChange              
            //   onCreateOption={(inputValue) => this.handleCreateOption('Selected Cell Types', inputValue)}
            //   options={options['Selected Cell Types']} // Set options to the fetched options
            />
          </div>



          {/* "Library Construction Method" (CreatableSelect) */}
          <div className="form-field">
            <label className="form-label">Library Construction Method:</label>
            <CreatableSelect
              name="Library Construction Method"
              value={formData['Library Construction Method']}
              onChange={handleChange}

            //   isClearable
            //   isSearchable
            //   isLoading={isLoading}
            //   onChange={(selectedOption) => this.handleSelectChange('Library Construction Method', selectedOption)} // Use handleSelectChange              
            //   onCreateOption={(inputValue) => this.handleCreateOption('Library Construction Method', inputValue)}
            //   options={options['Library Construction Method']} // Set options to the fetched options
              className="form-input"
            />
          </div>


          {/* "Nucleic Acid Source" (CreatableSelect) */}
          <div className="form-field">
            <label className="form-label">Nucleic Acid Source:</label>
            <CreatableSelect
              name="Nucleic Acid Source"
              value={formData['Nucleic Acid Source']}
              onChange={handleChange}

            //   isClearable
            //   isSearchable
            //   isLoading={isLoading}
            //   onChange={(selectedOption) => this.handleSelectChange('Nucleic Acid Source', selectedOption)} // Use handleSelectChange              
            //   onCreateOption={(inputValue) => this.handleCreateOption('Nucleic Acid Source', inputValue)}
            //   options={options['Nucleic Acid Source']} // Set options to the fetched options
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
                onChange={handleChange}

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
                  className="form-input"
                  onChange={handleChange}

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
            //   onChange={this.handleChange}
            onChange={handleChange}

              className="form-input"
            />
          </div>

          {/* "Disease Status (Specimen)" (CreatableSelect) */}
          <div className="form-field">
            <label className="form-label">Disease Status (Specimen):</label>
            <CreatableSelect
              name="Disease Status (Specimen)"
              value={formData['Disease Status (Specimen)']}
              onChange={handleChange}

            //   isClearable
            //   isSearchable
            //   isLoading={isLoading}
            //   onChange={(selectedOption) => this.handleSelectChange('Disease Status (Specimen)', selectedOption)} // Use handleSelectChange              
            //   onCreateOption={(inputValue) => this.handleCreateOption('Disease Status (Specimen)', inputValue)}
            //   options={options['Disease Status (Specimen)']} // Set options to the fetched options
            />
            {/* {errors['Disease Status (Specimen)'] && <div className="error-tooltip">{errors['Disease Status (Specimen)']}</div>} */}
          </div>


          {/* "Disease Status (Donor)" (CreatableSelect) */}
          <div className="form-field">
            <label className="form-label">Disease Status (Donor):</label>
            <CreatableSelect
              name="Disease Status (Donor)"
              value={formData['Disease Status (Donor)']}
              onChange={handleChange}

            //   isClearable
            //   isSearchable
            //   isLoading={isLoading}
            //   onChange={(selectedOption) => this.handleSelectChange('Disease Status (Donor)', selectedOption)} // Use handleSelectChange              
            //   onCreateOption={(inputValue) => this.handleCreateOption('Disease Status (Donor)', inputValue)}
            //   options={options['Disease Status (Donor)']} // Set options to the fetched options
            />
            {/* {errors['Disease Status (Donor)'] && <div className="error-tooltip">{errors['Disease Status (Donor)']}</div>} */}
          </div>

          {/* "Development Stage" (CreatableSelect) */}
          <div className="form-field">
            <label className="form-label">Development Stage:</label>
            <CreatableSelect
              name="Development Stage"
              value={formData['Development Stage']}
              onChange={handleChange}

            //   isClearable
            //   isSearchable
            //   isLoading={isLoading}
            //   onChange={(selectedOption) => this.handleSelectChange('Development Stage', selectedOption)} // Use handleSelectChange              
            //   onCreateOption={(inputValue) => this.handleCreateOption('Development Stage', inputValue)}
            //   options={options['Development Stage']} // Set options to the fetched options
              className="form-input"
            />
          </div>

          <div className="form-field">
            <label className="form-label">Donor Count:</label>
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
          <div className="form-field">
            <label className="form-label">Source:</label>
            <CreatableSelect
              name="Source"
              value={formData['Source']}
              onChange={handleChange}

            //   isClearable
            //   isSearchable
            //   isLoading={isLoading}
            //   onChange={(selectedOption) => this.handleSelectChange('Source', selectedOption)} // Use handleSelectChange              
            //   onCreateOption={(inputValue) => this.handleCreateOption('Source', inputValue)}
            //   options={options['Source']} // Set options to the fetched options
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
            //   onChange={this.handleChange}
            onChange={handleChange}

              placeholder="Enter ..."
              className="form-input"
            />
            {/* {errors['Source Key'] && <p className="error">{errors['Source Key']}</p>} */}
          </div>

          <div className="form-field">
            <label>Submission Date:</label>
            <input
              type="date"
              required
              name="Submission Date"
              value={FormData["Submission Date"]}
              onChange={handleChange}

            //   onChange={this.handleChange}
            />
            {/* {errors['Submission Date'] && <div className="error-tooltip">{errors['Submission Date']}</div>} */}
          </div>

          <div style={{ display: 'flex', justifyContent: 'center', marginTop: '20px' }}>
  <div style={{ marginRight: '10px' }}>
    <button type="submit" onClick={submitHandle}>Submit</button>
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
        </div>
        </div>
    );
};

export default ResultsTable;
