import React, { useEffect, useState } from 'react';
import { useLocation } from 'react-router-dom';
import axios from 'axios';
import { NODE_API_URL } from '../../../constants/declarations';
import { ScaleLoader } from 'react-spinners';
import AlertMessageComponent from '../../publishDatasets/components/alertMessageComponent';
import { Card, CardContent, Typography } from '@material-ui/core';
import { faAngleDown, faAngleRight } from '@fortawesome/free-solid-svg-icons';
import { FontAwesomeIcon } from '@fortawesome/react-fontawesome';
import { downloadFile, getFileNameFromURL, getCookie } from '../../../utils/utilFunctions';
import BenchmarksPlots from '../../publishDatasets/components/benchmarksPlots';
import RightRail from '../../RightNavigation/rightRail';
import DatasetDetailsTable from '../../Benchmarks/components/DatasetDetailsTable';


const DatasetInfoComponent = () => {
  const location = useLocation();
  const [datasetId, setDatasetId] = useState(null);
  const [datasetDetails, setDatasetDetails] = useState([]);
  const [loading, setLoading] = useState(true);
  const [message, setMessage] = useState('');
  const [hasMessage, setHasMessage] = useState(false);
  const [isError, setIsError] = useState(false);

  const [sectionsVisibility, setSectionsVisibility] = useState({
    metadataInfo: true,
    preprocessResults: true
  });

  const toggleSectionVisibility = (section) => {
    setSectionsVisibility((prevVisibility) => ({
      ...prevVisibility,
      [section]: !prevVisibility[section],
    }));
  };

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

  if (loading) {
    return (
      <div className="spinner-container">
        <ScaleLoader color="#36d7b7" loading={loading} />
      </div>
    );
  }

  return (
    <div className="single-dataset-details-page eighty-twenty-grid">
        <div className="main-content">
      {hasMessage && (
        <AlertMessageComponent message={message} setHasMessage={setHasMessage} setMessage={setMessage} isError={isError} />
      )}
      <h1>Details for dataset Id - {datasetId}</h1>
      {datasetDetails.length > 0 ? (
        datasetDetails.map((detail, index) => (
          <div key={index}>
            <div className='section'>
              <div className='section-heading' onClick={() => toggleSectionVisibility('metadataInfo')}>
                <h3>Dataset Metadata</h3>
                <span className="category-icon">
                  <FontAwesomeIcon
                    icon={sectionsVisibility.metadataInfo ? faAngleDown : faAngleRight}
                  />
                </span>
              </div>
              <div className='section-content' style={{ display: sectionsVisibility.metadataInfo ? 'block' : 'none' }}>
                <Card className="dataset-info-results">
                  <CardContent>
                    <DatasetDetailsTable 
                      datasetDetails={detail.datasetDetails} 
                      downloadFile={downloadFile} 
                      getFileNameFromURL={getFileNameFromURL} 
                    />
                  </CardContent>
                </Card>
              </div>
            </div>
            <div className='section'>
              <div className='section-heading' onClick={() => toggleSectionVisibility('preprocessResults')}>
                <h3>Benchmarks</h3>
                <span className="category-icon">
                  <FontAwesomeIcon
                    icon={sectionsVisibility.preprocessResults ? faAngleDown : faAngleRight}
                  />
                </span>
              </div>
              <div className='section-content' style={{ display: sectionsVisibility.preprocessResults ? 'block' : 'none' }}>
                <React.Fragment>
                  <Card className="benchmarks-results">
                    <CardContent>
                      <Typography variant="body2">PreProcess Results for {detail.datasetId}</Typography>
                      {/* <BenchmarksPlots
                        benchmarksPlot={detail.benchmarks_plot}
                        utilizationPlot={detail.utilization_plot}
                      /> */}
                    </CardContent>
                  </Card>
                </React.Fragment>
              </div>
            </div>
          </div>
        ))
      ) : (
        <p>No data found</p>
      )}
      </div>

    <div>
    {(getCookie('jwtToken') != undefined || getCookie('jwtToken') != '') && (<div className="right-rail">
        <RightRail />
    </div>)}
      </div>
    </div>
  );
};

export default DatasetInfoComponent;
