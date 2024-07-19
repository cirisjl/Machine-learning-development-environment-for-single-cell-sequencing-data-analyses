import React, { useEffect, useState } from 'react';
import { useLocation } from 'react-router-dom';
import axios from 'axios';
import { NODE_API_URL } from '../../../constants/declarations';
import { ScaleLoader } from 'react-spinners';
import AlertMessageComponent from '../../publishDatasets/components/alertMessageComponent';
import { Card, CardContent, Typography} from '@material-ui/core';
import { faAngleDown, faAngleRight } from '@fortawesome/free-solid-svg-icons';
import { FontAwesomeIcon } from '@fortawesome/react-fontawesome';
import { downloadFile, getFileNameFromURL } from '../../../utils/utilFunctions';
import BenchmarksPlots from '../../publishDatasets/components/benchmarksPlots';

const BenchmarksViewDetailsComponent = () => {
    const location = useLocation();
    const [benchmarksId, setBenchmarksId] = useState(null);
    const [benchmarksDetails, setBenchmarksDetails] = useState([]);
    const [loading, setLoading] = useState(true);
    const [message, setMessage] = useState('');
    const [hasMessage, setHasMessage] = useState(false);
    const [isError, setIsError] = useState(false);

    const [sectionsVisibility, setSectionsVisibility] = useState({
        taskInfo: true,
        benchmarks: true,
    });

    const toggleSectionVisibility = (section) => {
        setSectionsVisibility((prevVisibility) => ({
          ...prevVisibility,
          [section]: !prevVisibility[section],
        }));
    };

    useEffect(() => {
        const queryParams = new URLSearchParams(location.search);
        const id = queryParams.get('benchmarksId');
        setBenchmarksId(id);

        if (id) {
            // Fetch the dataset based on benchmarksId using axios
            axios.post(`${NODE_API_URL}/getBenchmarksResults`, { benchmarksId: id })
                .then(response => {
                    setBenchmarksDetails(response.data);
                    setMessage(`Successfully fetched details for the benchmarks ID - ${id}.`);
                    setHasMessage(true);
                    setIsError(false);
                    setLoading(false);
                })
                .catch(error => {
                    console.error(`Error fetching benchmarks details for benchmark ID - ${id}:`, error);
                    setMessage(`Error fetching benchmarks details for benchmark ID - ${id}.`);
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
        <div className="single-benchmark-details-page">
            {hasMessage && (
                <AlertMessageComponent message={message} setHasMessage={setHasMessage} setMessage={setMessage} isError={isError} />
            )}
            <h1>Details for BenchmarksId - {benchmarksId}</h1>
            {benchmarksDetails.length > 0 ? (
                benchmarksDetails.map((detail, index) => (
                    <div key={index}>
                        <div className='section'>
                            <div className='section-heading' onClick={() => toggleSectionVisibility('taskInfo')}>
                                <h3>Task Info</h3>
                                <span className="category-icon">
                                    <FontAwesomeIcon
                                        icon={sectionsVisibility.taskInfo ? faAngleDown : faAngleRight}
                                    />
                                </span>
                            </div>
                            <div className='section-content' style={{ display: sectionsVisibility.taskInfo ? 'block' : 'none' }}>
                                <Card className="benchmarks-results">
                                    <CardContent>
                                        <div key={detail.benchmarksId} className="dataset-details">
                                            <Typography variant="body2">
                                            <ul>
                                                    <li><strong>Dataset ID:</strong> {detail.datasetId}</li>
                                                    <li><strong>Task Type:</strong> {detail.task_type}</li>
                                                    <li><strong>Archive Path:</strong> 
                                                        <a 
                                                            download 
                                                            onClick={() => { downloadFile(detail.archive_path) }} 
                                                            style={{ 
                                                                marginLeft: '10px', 
                                                                textAlign: 'center', 
                                                                cursor: 'pointer', 
                                                                textDecoration: 'underline', 
                                                                color: 'blue' 
                                                            }}>
                                                            {getFileNameFromURL(detail.archive_path) || 'Not available'}
                                                        </a>
                                                    </li>
                                                    <li><strong>Train Path:</strong> 
                                                        <a 
                                                            download 
                                                            onClick={() => { downloadFile(detail.train_path) }} 
                                                            style={{ 
                                                                marginLeft: '10px', 
                                                                textAlign: 'center', 
                                                                cursor: 'pointer', 
                                                                textDecoration: 'underline', 
                                                                color: 'blue' 
                                                            }}>
                                                            {getFileNameFromURL(detail.train_path) || 'Not available'}
                                                        </a>
                                                    </li>
                                                    <li><strong>Test Path:</strong> 
                                                        <a 
                                                            download 
                                                            onClick={() => { downloadFile(detail.test_path) }} 
                                                            style={{ 
                                                                marginLeft: '10px', 
                                                                textAlign: 'center', 
                                                                cursor: 'pointer', 
                                                                textDecoration: 'underline', 
                                                                color: 'blue' 
                                                            }}>
                                                            {getFileNameFromURL(detail.test_path) || 'Not available'}
                                                        </a>
                                                    </li>
                                                    <li><strong>Validation Path:</strong> 
                                                        <a 
                                                            download 
                                                            onClick={() => { downloadFile(detail.validation_path) }} 
                                                            style={{ 
                                                                marginLeft: '10px', 
                                                                textAlign: 'center', 
                                                                cursor: 'pointer', 
                                                                textDecoration: 'underline', 
                                                                color: 'blue' 
                                                            }}>
                                                            {getFileNameFromURL(detail.validation_path) || 'Not available'}
                                                        </a>
                                                    </li>
                                                </ul>
                                            </Typography>
                                        </div>
                                    </CardContent>
                                </Card>
                            </div>
                        </div>

                        <div className='section'>
                            <div className='section-heading' onClick={() => toggleSectionVisibility('benchmarks')}>
                                <h3>Benchmarks</h3>
                                <span className="category-icon">
                                    <FontAwesomeIcon
                                        icon={sectionsVisibility.benchmarks ? faAngleDown : faAngleRight}
                                    />
                                </span>
                            </div>
                            <div className='section-content' style={{ display: sectionsVisibility.benchmarks ? 'block' : 'none' }}>
                                <React.Fragment>
                                    <Card className="benchmarks-results">
                                        <CardContent>
                                            <Typography variant="body2">Benchmark Results for {detail.datasetId}</Typography>
                                            <BenchmarksPlots
                                                benchmarksPlot={detail.benchmarks_plot}
                                                utilizationPlot={detail.utilization_plot}
                                            />
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
    );
};

export default BenchmarksViewDetailsComponent;
