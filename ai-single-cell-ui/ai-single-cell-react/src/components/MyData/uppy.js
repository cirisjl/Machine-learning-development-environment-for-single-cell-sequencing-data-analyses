import React, { useEffect, useState } from 'react'
import Uppy from '@uppy/core'
import XHRUpload from '@uppy/xhr-upload'
import { Dashboard } from '@uppy/react'
import GoogleDrive from '@uppy/google-drive'
import OneDrive from '@uppy/onedrive'
import Dropbox from '@uppy/dropbox'
import "@uppy/dashboard/dist/style.css"
import "@uppy/core/dist/style.css"
import "@uppy/progress-bar/dist/style.css"
import "@uppy/status-bar/dist/style.css"
import "@uppy/drag-drop/dist/style.css"

const SERVER_URL = "http://" + process.env.REACT_APP_HOST_URL + ":3001";
export default function UppyUploader(props) {
    const { isUppyModalOpen, setIsUppyModalOpen, pwd, authToken, freeSpace } = props;
    const [dimensions, setDimensions] = useState({ width: 0, height: 0 });

    useEffect(() => {
        setDimensions({
            width: window.innerWidth,
            height: window.innerHeight,
        });
    }, [window.innerWidth, window.innerHeight]);


    useEffect(() => {
        function handleResize() {
            setDimensions({
                width: window.innerWidth,
                height: window.innerHeight,
            });
        }

        window.addEventListener('resize', handleResize);

        // Cleanup function that removes the event listener
        return () => window.removeEventListener('resize', handleResize);
    }, []);

    const uppy = new Uppy({
        id: 'fileUploader',
        autoProceed: false,
        allowMultipleUploads: true,
        restrictions: {
            maxFileSize: freeSpace * 1024 * 1024 * 1024,
            maxNumberOfFiles: 5,
            maxTotalFileSize: freeSpace * 1024 * 1024 * 1024,
        },
        debug: true,
    });
    uppy.use(GoogleDrive, {
        companionUrl: `http://${process.env.REACT_APP_HOST_URL}:3020`,
    });
    uppy.use(OneDrive, {
        companionUrl: `http://${process.env.REACT_APP_HOST_URL}:3020`,
    });
    uppy.use(Dropbox, {
        companionUrl: `http://${process.env.REACT_APP_HOST_URL}:3020`,
    });
    uppy.use(XHRUpload, {
        endpoint: `${SERVER_URL}/upload?uploadDir=${pwd}&authToken=${authToken}`,
        formData: true,
        fieldName: 'files'
    });
    uppy.on('file-added', (file) => {
        console.log('Upload destination:', file.meta.destination);
        uppy.upload();
    });
    if (isUppyModalOpen)
        return (<div className="uppy-modal">
            <Dashboard uppy={uppy} plugins={['GoogleDrive', 'OneDrive', 'Dropbox']} />
            <button style={{
                top: `${dimensions.height * 0.5 + 245}px`,
                left: `${dimensions.width * 0.5 + 330}px`,
                position: "absolute",
                transform: "translate(-50%, -50%)",
                padding: "5px 5px",
                cursor: "pointer",
                border: "1px solid black",
                borderRadius: "3px"
            }}
                onClick={() => { setIsUppyModalOpen(!isUppyModalOpen) }}
            >Close
            </button>
        </div>
        )
}


