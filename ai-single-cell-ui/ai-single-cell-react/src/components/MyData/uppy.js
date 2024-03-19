import React, { useEffect, useState } from 'react'
import Uppy from '@uppy/core'
import XHRUpload from '@uppy/xhr-upload'
import { Dashboard } from '@uppy/react'
import GoogleDrive from '@uppy/google-drive'
import OneDrive from '@uppy/onedrive'
import Dropbox from '@uppy/dropbox'
import Url from '@uppy/url';
import "@uppy/dashboard/dist/style.css"
import "@uppy/core/dist/style.css"
import "@uppy/progress-bar/dist/style.css"
import "@uppy/status-bar/dist/style.css"
import "@uppy/drag-drop/dist/style.css"
import { setHashStore, getHashStore, calcMD5Hash } from '../../utils/utilFunctions'

import axios from 'axios';

const SERVER_URL = "http://" + process.env.REACT_APP_HOST_URL + ":3001";

const checkIfFileExists = async (file) => {
    const hashes = getHashStore().datasetHashes;
    console.log(hashes)

    if (hashes) {
        for (let item of hashes) {
            if (item.hasOwnProperty('hash') && item.hash === await calcMD5Hash(file.data)) {
                return true;
            }
        }
        return false;
    }
    return true;
}

export default function UppyUploader(props) {

    const { isUppyModalOpen, setIsUppyModalOpen, pwd, authToken, freeSpace, publicDatasetFlag, toPublishDataset, setFileError, setTaskData, context } = props;
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
    uppy.use(Url, {
        companionUrl: `http://${process.env.REACT_APP_HOST_URL}:3020`,
    });
    uppy.use(XHRUpload, {
        endpoint: `${SERVER_URL}/upload?uploadDir=${pwd}&authToken=${authToken}&publicDatasetFlag=${publicDatasetFlag}`,
        formData: true,
        fieldName: 'files'
    });

    uppy.on('file-added', async (file) => {
        let hash = await calcMD5Hash(file.data);
        console.log("hash")
        console.log(hash)

        // Check if the hash for the file already exists in the system
        let res = await axios.get(`${SERVER_URL}/mongoDB/api/file-exists?hash=${hash.hashResult}`)
        console.log(res.data)

        // Add the hash to global store so it can be added as dataset metadata
        setHashStore({ ...getHashStore(), datasetHashes: [...(getHashStore().datasetHashes || []), { file: file.name, hash: hash.hashResult, exists: res.data.exists }] })

        // Remove the file so it can't be uploaded
        if (res.data.exists === true) {
            console.log(`Removed file ${file.name}`)
            uppy.info(`Removed file ${file.name} because it already exists.`)
            uppy.removeFile(file.id)
        }
    });

    uppy.on('upload-success', async (file, response) => {
        if (toPublishDataset) {
            // Access the filename of the successfully uploaded file
            setFileError('');

            const hash = await calcMD5Hash(file.data)
            console.log("hash from upload-success inside if")
            console.log(hash)

            setTaskData((prevTaskData) => ({
                ...prevTaskData,
                upload: {
                    ...prevTaskData.upload,
                    files: [...(prevTaskData.upload.files || []), file.name],
                },
            }));
            console.log('Successfully uploaded file name:', file.name);
        }

        // Generate and save file hash on upload success
        const hash = await calcMD5Hash(file.data)

        // Add file info to file collection for dedup
        let res = await axios.get(`${SERVER_URL}/mongoDB/api/file-exists?hash=${hash.hashResult}`)
        console.log(res.data)

        if (res.data.exists === false) {
            const payload = {
                name: file.name,
                size: file.size,
                type: file.type,
                hash: hash.hashResult,
                uploadedBy: authToken,
                createdAt: new Date().toUTCString()
            }
            await axios.post(`${SERVER_URL}/mongoDB/api/add-file?authToken=${authToken}`, payload).then(res => { console.log(res) })
        }
    });

    if (isUppyModalOpen && !toPublishDataset)
        return (<div className="uppy-modal">
            <Dashboard uppy={uppy} plugins={['GoogleDrive', 'OneDrive', 'Dropbox', 'Url']} />
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

    if (toPublishDataset) {
        return (<div className='uppy-comp'>
            <Dashboard uppy={uppy} plugins={['GoogleDrive', 'OneDrive', 'Dropbox', 'Url']} />
        </div>)
    }
}


