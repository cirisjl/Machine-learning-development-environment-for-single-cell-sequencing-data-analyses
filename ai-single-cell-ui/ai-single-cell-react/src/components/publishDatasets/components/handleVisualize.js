import React, { useState, useEffect } from 'react';
import { useParams } from 'react-router-dom';
import { SERVER_URL } from '../../../constants/declarations';
import axios from 'axios';
import '../publishDatasets.css';
import ReactPlotly from './reactPlotly';



function HandleVisualize() {
  const { id } = useParams();
  const [datasetId, setDatasetId] = useState("");
  const [data, setData] = useState(null);
  const [visualizeDataset,setVisualizeDataset]=useState({});
  const [ppresult,setppResult]=useState({});

  const [activeTab, setActiveTab] = useState(0);
  const [selectedPlot, setSelectedPlot] = useState('2d'); // State to manage selected plot option

  const handlePlotChange = (event) => {
    setSelectedPlot(event.target.value);
    console.log('check',ppresult);

  };

  const handleTabClick = (tabIndex) => {
    setActiveTab(tabIndex);
  };
  
  useEffect(() => {
    console.log(" id:", id);
    setDatasetId(id);
  }, [id]);

 

  useEffect(() => {
    // Make the API call here
    const fetchData = async () => {
    
      axios.get(`${SERVER_URL}/mongoDB/api/fetchDataforVisualize/${id}`)
      .then(response => {
        console.log('Response Received',response.data);
        setVisualizeDataset(response.data)
        
      })
      .catch(error => {
        console.log('Error while retrieving data',error);
        
      });
      
    
    }

    fetchData();
  }, [id]);

  useEffect(() => {
    if(visualizeDataset){
    axios.get(`${SERVER_URL}/mongoDB/api/fetchGraphData/${visualizeDataset.process_id}`)
      .then(response => {
        console.log('Response Received for pp dataset',response.data);
        setppResult(response.data)

      })
      .catch(error => {
        console.log('Error while retrieving data',error);
      });
    }
    
  }, [visualizeDataset]);

  return (


<div className="dataset-container">
  <div className="dataset-dialog-container">
    
    <div className="dataset-dialog">
      <div className="dataset-tabs">
        <button className={activeTab === 0 ? 'active-tab' : ''} onClick={() => handleTabClick(0)}>Summary</button>
        <button className={activeTab === 1 ? 'active-tab' : ''} onClick={() => handleTabClick(1)}>Explore</button>
        <button className={activeTab === 2 ? 'active-tab' : ''} onClick={() => handleTabClick(2)}>Download</button>
      </div>
      <div className="dataset-tab-content">
        {activeTab === 0 && Object.keys(visualizeDataset).length !== 0 && (
          <div className="dataset-tab-panel">
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
          </div>
        )}
        {/* {activeTab === 1 &&
          <div className="dataset-tab-panel">
   <h2>Scatter Plot</h2>
           <ReactPlotly plot_data={ppresult.scatter_plot} />

           <h2>UMAP Plot</h2>
           <ReactPlotly plot_data={ppresult.umap_plot} />


           <h2>Violin Plot</h2>
           <ReactPlotly plot_data={ppresult.violin_plot} />

           <h2>3d plot</h2>
           {ppresult.umap_plot_3d ? (
        <ReactPlotly plot_data={ppresult.umap_plot_3d} />
      ) : (
        <div>No 3D plot data available</div>
      )}
          </div>
        } */}

{activeTab === 1 &&
              <div className="dataset-tab-panel">
                <div id="plot-container"> 
                  <label htmlFor="plot-select">Select plot type: </label>
                  <select id="plot-select" value={selectedPlot} onChange={handlePlotChange}>
                    <option value="2d">2D Plots</option>
                    <option value="3d">3D Plot</option>
                  </select>
                </div>
                {selectedPlot === '2d' ? (
                  <>
                  <h2>UMAP Plot</h2>
                  {ppresult.umap_plot ? (
        <ReactPlotly plot_data={ppresult.umap_plot} />
          ) : (
          <div>No UMAP plot data available</div>
           )}
                    <h2>Scatter Plot</h2>
                    {ppresult.scatter_plot ? (
      <ReactPlotly plot_data={ppresult.scatter_plot} />
    ) : (
      <div>No scatter plot data available</div>
    )}
                    <h2>Violin Plot</h2>

                    {ppresult.violin_plot ? (
      <ReactPlotly plot_data={ppresult.violin_plot} />
    ) : (
      <div>No violin plot data available</div>
    )}
                  </>
                ) : (
                  <>
                    <h2>3D Plot</h2>
                    {ppresult.umap_plot_3d ? (
                      <ReactPlotly plot_data={ppresult.umap_plot_3d} />
                    ) : (
                      <div>No 3D plot data available</div>
                    )}
                  </>
                )}
              </div>
            }
        {activeTab === 2 &&
          <div className="dataset-tab-panel">
            Content for Tab 3
          </div>
        }
      </div>
    </div>
    {activeTab === 0 && Object.keys(visualizeDataset).length !== 0 && (
      <div className="side-panel">
        <h2>Species</h2>
          <p>{visualizeDataset['Species'].value}</p>
           <hr/>
           <h2>Sample Type</h2>
           <p>{visualizeDataset['Sample Type'].value}</p>
           <hr/>
           <h2>Organ Part</h2>
           <p>{visualizeDataset['Organ Part'].label}</p>
           <hr/>
           <h2>Model Organ</h2>
           <p>{visualizeDataset['Model Organ'].label}</p>
           <hr/>
           <h2>Selected Cell Types</h2>
          <p>{visualizeDataset['Selected Cell Types'].label}</p>
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
           

      </div>
    )}
  </div>
</div>


  );
}

export default HandleVisualize;