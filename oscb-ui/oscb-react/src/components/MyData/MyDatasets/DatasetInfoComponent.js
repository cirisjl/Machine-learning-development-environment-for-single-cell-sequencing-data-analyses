import React, { useEffect, useState } from 'react';
import { useLocation } from 'react-router-dom';
import axios from 'axios';
import { NODE_API_URL } from '../../../constants/declarations';
import { ScaleLoader } from 'react-spinners';
import AlertMessageComponent from '../../publishDatasets/components/alertMessageComponent';
import {Card, CardContent, Typography ,List, ListItem, ListItemText } from '@mui/material';
import { faAngleDown, faAngleRight } from '@fortawesome/free-solid-svg-icons';
import { FontAwesomeIcon } from '@fortawesome/react-fontawesome';
import { downloadFile, getFileNameFromURL, getCookie } from '../../../utils/utilFunctions';
import RightRail from '../../RightNavigation/rightRail';
import DatasetDetailsTable from '../../Benchmarks/components/DatasetDetailsTable';

import { styled } from '@mui/material/styles';
import ArrowForwardIosSharpIcon from '@mui/icons-material/ArrowForwardIosSharp';
import MuiAccordion from '@mui/material/Accordion';
import MuiAccordionSummary from '@mui/material/AccordionSummary';
import MuiAccordionDetails from '@mui/material/AccordionDetails';

const DatasetInfoComponent = () => {
  const location = useLocation();
  const [datasetId, setDatasetId] = useState(null);
  const [datasetDetails, setDatasetDetails] = useState([]);
  const [loading, setLoading] = useState(true);
  const [message, setMessage] = useState('');
  const [hasMessage, setHasMessage] = useState(false);
  const [isError, setIsError] = useState(false);
  const [expanded, setExpanded] = React.useState(0);

  const Accordion = styled((props) => (
    <MuiAccordion disableGutters elevation={0} square {...props} />
  ))(({ theme }) => ({
    border: `1px solid ${theme.palette.divider}`,
    '&:not(:last-child)': {
      borderBottom: 0,
    },
    '&::before': {
      display: 'none',
    },
  }));
  
  const AccordionSummary = styled((props) => (
    <MuiAccordionSummary
      expandIcon={<ArrowForwardIosSharpIcon sx={{ fontSize: '0.9rem' }} />}
      {...props}
    />
  ))(({ theme }) => ({
    backgroundColor: 'rgba(0, 0, 0, .03)',
    flexDirection: 'row-reverse',
    '& .MuiAccordionSummary-expandIconWrapper.Mui-expanded': {
      transform: 'rotate(90deg)',
    },
    '& .MuiAccordionSummary-content': {
      marginLeft: theme.spacing(1),
    },
    ...theme.applyStyles('dark', {
      backgroundColor: 'rgba(255, 255, 255, .05)',
    }),
  }));
  
  const AccordionDetails = styled(MuiAccordionDetails)(({ theme }) => ({
    padding: theme.spacing(2),
    borderTop: '1px solid rgba(0, 0, 0, .125)',
  }));

  const handleChange = (panel) => (event, newExpanded) => {
    setExpanded(newExpanded ? panel : false);
  };

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
                <h3>Pre Process Results</h3>
                <span className="category-icon">
                  <FontAwesomeIcon
                    icon={sectionsVisibility.preprocessResults ? faAngleDown : faAngleRight}
                  />
                </span>
              </div>
              <div className='section-content' style={{ display: sectionsVisibility.preprocessResults ? 'block' : 'none' }}>
                <React.Fragment>
                  <Card className="dataset-info-results">
                    <CardContent>
                      {/* <Typography variant="body2">PreProcess Results for {detail.datasetId}</Typography> */}
                        <div>
                            {detail.preProcessResults.length > 0 ? (
                                detail.preProcessResults.map((preProcessResult, index) => (
                                    <Accordion expanded={expanded === index} onChange={handleChange(index)}>
                                        <AccordionSummary>
                                            <Typography>{preProcessResult.description}</Typography>
                                        </AccordionSummary>
                                        <AccordionDetails>
                                            <Descriptions title="Pre Process Result" bordered column={1}>
                                                <Descriptions.Item label="Description">{preProcessResult.description}</Descriptions.Item>
                                                <Descriptions.Item label="Stage">{preProcessResult.stage}</Descriptions.Item>
                                                <Descriptions.Item label="Process">{preProcessResult.process}</Descriptions.Item>
                                                <Descriptions.Item label="Method">{preProcessResult.method}</Descriptions.Item>
                                                <Descriptions.Item label="Parameters">
                                                <Card sx={{ maxWidth: 600, margin: '20px auto' }}>
                                                    <CardContent>
                                                        <List>
                                                        {Object.entries(preProcessResult.parameters).map(([key, value]) => (
                                                            <ListItem key={key} divider>
                                                            <ListItemText
                                                                primary={key}
                                                                secondary={String(value)}
                                                                primaryTypographyProps={{ variant: 'body1', fontWeight: 'bold' }}
                                                                secondaryTypographyProps={{ variant: 'body2', color: 'text.secondary' }}
                                                            />
                                                            </ListItem>
                                                        ))}
                                                        </List>
                                                    </CardContent>
                                                </Card>
                                                </Descriptions.Item>
                                                <Descriptions.Item label="Files"></Descriptions.Item>
                                            </Descriptions>
                                            <div>
                                                <p>Plots: </p>
                                            </div>
                                        </AccordionDetails>
                                    </Accordion>
                                ))
                            )  : (
                                <p>No Pre Process Results to show for the dataset.</p>
                            )}              
                        </div>
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