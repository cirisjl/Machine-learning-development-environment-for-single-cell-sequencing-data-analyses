const express = require('express');
const fs = require('fs-extra');
const mysql = require('mysql2');
const bcrypt = require('bcrypt');
const jwt = require('jsonwebtoken');
const bodyParser = require('body-parser');
const cors = require('cors');
const cookieParser = require('cookie-parser');
const path = require('path')
const mime = require('mime');
const archiver = require('archiver');
const util = require('util');
const stat = util.promisify(fs.stat);
const multer = require("multer");
const hostIp = process.env.SSH_CONNECTION.split(' ')[2];
require('dotenv').config();

const mongoDBConfig = JSON.parse(fs.readFileSync('./configs/mongoDB.json'));// Import the MongoDB connection configuration
const { mongoUrl, dbName, optionsCollectionName, datasetCollection, taskBuilderCollectionName, userDatasetsCollection } = mongoDBConfig;
const { MongoClient, ObjectId } = require('mongodb');

// const Option = require('../models/Option');
// // Import the database configuration
// require('./config/mongoDBClient');

// Increase the limit for the request body size to 25MB

console.log('HOSTURL: ' + process.env.HOST_URL);
const app = express();
app.use(cors({
    origin: [`http://${process.env.HOST_URL}:3000`, `http://${hostIp}:3000`, 'http://node-0.jiangl0-160204.biomizzou-pg0.clemson.cloudlab.us:3000', 'http://node-0.jiangl0-161295.biomizzou-pg0.clemson.cloudlab.us:3000', 'http://node-0.ai-single-cell.biomizzou-pg0.clemson.cloudlab.us:3000', 'http://node-1.c220g2-sampath.biomizzou-pg0.wisc.cloudlab.us:3000'],
    credentials: true
}));
app.use(bodyParser.json({ limit: '25mb' }));
app.use(cookieParser());
app.use((err, req, res, next) => {
    console.error(err.stack);
    res.status(500).send('Unhandled Error in Node Application');
});

const dbConfig = JSON.parse(fs.readFileSync('./configs/dbconfigs.json'));
const storageConfig = JSON.parse(fs.readFileSync('./configs/storageConfig.json'));
const { storageDir, storageAllowance, intermediateStorage, publicStorage } = storageConfig;

// Create a connection pool to handle multiple connections to the database
const pool = mysql.createPool({
    host: dbConfig.host,
    user: dbConfig.user,
    password: dbConfig.password,
    database: dbConfig.database,
    connectionLimit: dbConfig.connectionLimit
});

// Middleware function to verify JWT token
function verifyToken(req, res, next) {
    const bearerHeader = req.headers['authorization'];

    if (typeof bearerHeader !== 'undefined') {
        const bearerToken = bearerHeader.split(' ')[1];
        req.token = bearerToken;
        next();
    } else {
        res.sendStatus(403);
    }
}

function getUserFromToken(token) {
    if (typeof token !== 'string') {
        return 'Unauthorized';
    }

    try {
        const decoded = jwt.verify(token, process.env.JWT_TOKEN_SECRET);
        if (!decoded.username) {
            return 'Unauthorized';
        }

        return decoded.username;
    } catch (err) {
        console.log('Session Expired. Please login again ' + err)
        return 'Unauthorized';
    }
}
const createDirectoryIfNotExists = async (dirPath) => {
    try {
        await fs.mkdir(dirPath, { recursive: true });
        console.log(`Directory "${dirPath}" created successfully.`);
    } catch (err) {
        if (err.code !== 'EEXIST') {
            console.error('Error creating the directory:', err);
        }
    }
};

const createUniqueFolder = (destinationDir, folderName, index = 1) => {
    const targetFolderName = index === 1 ? path.join(destinationDir, folderName) : path.join(destinationDir, `${folderName}(${index})`);
    const targetFolderPath = path.join(__dirname, targetFolderName);

    if (!fs.existsSync(targetFolderPath)) {
        try {
            fs.mkdirSync(targetFolderPath);
            console.log(`Directory "${targetFolderName}" created successfully.`);
            return targetFolderName;
        } catch (err) {
            console.error('Error creating the directory:', err);
            return null;
        }
    } else {
        return createUniqueFolder(destinationDir, folderName, index + 1); // Try with the next index
    }
};

// Function to copy files from source directory to destination directory
const copyFiles = async (sourceDir, destinationDir, dirName, files, fromPublic) => {
    try {


        for (let file of files) {
            const sourceFilePath = path.join(sourceDir, file);
            let destinationFilePath = "";
            if (fromPublic) {
                file = file.replace(/^\/?publicDatasets\//, '/'); // Remove "PUBLIC_DATASETS" from the start
                destinationFilePath = path.join(destinationDir, file);
            } else {
                destinationFilePath = path.join(destinationDir, file);
            }

            const sourceFileDir = path.dirname(sourceFilePath);
            const destinationFileDir = path.dirname(destinationFilePath);

            // Ensure the destination directory exists before copying files
            await createDirectoryIfNotExists(destinationFileDir);

            // Perform the actual file copy
            await fs.copyFile(sourceFilePath, destinationFilePath);
        }
    } catch (error) {
        console.error('Error copying files:', error);
    }
};

app.post('/api/copyFiles', async (req, res) => {

    const { selectedFiles, userId } = req.body;

    try {
        let filesFromPublic = false;
        let dirName = ""

        // Logic to Copy files from public storage to user private storage if it is a public Dataset.
        for (const file of selectedFiles) {
            if (file.startsWith("publicDataset") || file.startsWith("/publicDatasets")) {
                filesFromPublic = true;
                break;
            }
        }

        if (filesFromPublic) {

            if (selectedFiles.length > 0) {
                dirName = path.dirname(selectedFiles[0])
            }

            let userPrivateStorageDir = storageDir + userId // Change this to the user's private storage path

            // Copy files from public dataset directory to user's private storage
            copyFiles("/usr/src/app/storage/", userPrivateStorageDir, dirName, selectedFiles, filesFromPublic);
            res.json({ status: 200, message: 'Files copied successfully' });
        }
    } catch (error) {
        res.json({ status: 500, message: 'Error while copying files from source to destination' });
    }
});

// Refresh token endpoint
app.get('/api/refresh-token', verifyToken, (req, res) => {
    jwt.verify(req.token, process.env.JWT_TOKEN_SECRET, (err, authData) => {
        if (err) {
            res.sendStatus(403);
        } else {
            if (authData.username !== null && authData.username !== undefined) {
                const newToken = jwt.sign({ username: authData.username }, process.env.JWT_TOKEN_SECRET, { expiresIn: '2m' });
                res.cookie('jwtToken', newToken, { maxAge: 2 * 60 * 1000, path: "/" });
                res.json({ status: 200, message: 'Token refreshed', token: newToken });
            }
        }
    })
});

// Route to handle user signup
app.post('/api/signup', (req, res) => {
    const { username, email, password } = req.body;

    // Hash the password using bcrypt
    bcrypt.hash(password, 10, (err, hash) => {
        if (err) {
            console.error(err);
            res.json({ status: 500, message: 'Internal Server Error' });
            return;
        }

        // Insert the user into the database
        pool.query('INSERT INTO users (username, email, password_hash) VALUES (?, ?, ?)', [username, email, hash], (err) => {
            if (err) {
                console.error(err);
                res.json({ status: 500, message: 'An Account is already present with the given email or username. Please try to login or create a new account using different email.' });
                return;
            }

            try {
                if (!err) {
                    if (!fs.existsSync(storageDir + username))
                        fs.promises.mkdir(storageDir + username);

                    // Create JWT token and send it back to the client
                    const jwtToken = jwt.sign({ username }, process.env.JWT_TOKEN_SECRET, { expiresIn: '2m' });

                    // the cookie will be set with the name "jwtToken" and the value of the token
                    // the "httpOnly" and "secure" options help prevent XSS and cookie theft
                    // the "secure" option is only set if the app is running in production mode
                    // set the cookie with the JWT token on the response object
                    res.cookie("jwtToken", jwtToken, {
                        //httpOnly: true,
                        maxAge: 2 * 60 * 1000,
                        path: "/"
                        //secure: process.env.NODE_ENV === "production",
                    });
                    res.json({ status: 200, message: 'User account created successfully' });
                }
            }
            catch (error) {
                res.json({ status: 500, message: 'Error occured while creating a storage directory for the user' });
            }

        });
    });
});

// Route to handle user login
app.post('/api/login', (req, res) => {
    const { username, password } = req.body;

    pool.query('SELECT * FROM users WHERE username = ?', [username], (err, results) => {
        if (err) {
            console.error(err);
            res.json({ status: 500, message: 'Internal Server Error' });
            return;
        }

        if (results.length === 0) {
            res.json({ status: 401, message: 'Invalid credentials' });
            return;
        }

        const user = results[0];

        bcrypt.compare(password, user.password_hash, (err, isMatch) => {
            if (err) {
                console.error(err);
                res.json({ status: 500, message: 'Internal Server Error' });
                return;
            }

            if (!isMatch) {
                res.json({ status: 401, message: 'Invalid credentials' });
                return;
            }

            // Create JWT token and send it back to the client
            const jwtToken = jwt.sign({ username }, process.env.JWT_TOKEN_SECRET, { expiresIn: '2m' });

            // the cookie will be set with the name "jwtToken" and the value of the token
            // the "httpOnly" and "secure" options help prevent XSS and cookie theft
            // the "secure" option is only set if the app is running in production mode
            // set the cookie with the JWT token on the response object
            res.cookie("jwtToken", jwtToken, {
                //httpOnly: true,
                maxAge: 2 * 60 * 1000,
                path: "/"
                //secure: process.env.NODE_ENV === "production",
            });

            res.json({ status: 200, message: 'Logged in successfully', jwtToken });
        });
    });
});

// Route to handle protected resource
app.get('/protected', verifyToken, (req, res) => {
    jwt.verify(req.token, process.env.JWT_TOKEN_SECRET, (err, authData) => {
        if (err) {
            res.sendStatus(403);
        } else {
            if (authData.username !== null && authData.username !== undefined) {
                pool.query('SELECT isAdmin FROM users WHERE username = ?', authData.username, (err, results) => {
                    if (err) {
                        console.error(err);
                        res.json({ message: 'Internal Server Error' });
                        return;
                    }

                    if (results.length === 0) {
                        res.json({ message: 'Invalid credentials' });
                        return;
                    }

                    const adminFlag = results[0].isAdmin;

                    authData.isAdmin = (adminFlag == 1) ? true : false;

                    res.json({ message: 'You have access to the protected resource', authData });
                });
            }
        }
    });
});

app.post('/createDataset', async (req, res) => {

    const { title, n_cells, reference, summary, authToken, files, makeItpublic } = req.body;
    const username = getUserFromToken(authToken);

    let filesFromPublic = false;

    // Logic to Copy files from public storage to user private storage if it is a public Dataset.
    for (const file of files) {

        if (file.startsWith("publicDataset") || file.startsWith("/publicDatasets")) {
            filesFromPublic = true;
            break;
        }
    }
    if (filesFromPublic) {
        let dirName = ""

        if (files.length > 0) {
            dirName = path.dirname(files[0])
        }

        let userPrivateStorageDir = storageDir + username // Change this to the user's private storage path

        // Copy files from user's private storage to public dataset directory
        await copyFiles("/usr/src/app/storage/", userPrivateStorageDir, dirName, files, filesFromPublic);
    }

    pool.getConnection(function (err, connection) {
        if (err) {
            console.error('Error getting DB connection:', err);
            return res.status(500).json({ message: 'Database connection error' });
        }

        connection.beginTransaction(function (err) {
            if (err) {
                console.error('Error starting transaction:', err);
                connection.release();
                return res.status(500).json({ message: 'Transaction error' });
            }

            // Run SELECT command
            connection.query('SELECT user_id FROM users WHERE username = ? LIMIT 1', [username], function (err, userRows) {
                if (err) {
                    connection.rollback(function () {
                        connection.release();
                    });
                }

                const userId = userRows[0].user_id;

                if (!userId) {
                    res.status(400).send('User not found');
                    connection.rollback(function () {
                        connection.release();
                    });
                    return;
                }

                connection.query('INSERT INTO dataset (title, n_cells, reference, summary, user_id) VALUES (?, ?, ?, ?, ?)', [title, n_cells, reference, summary, userId], function (err, datasetResult) {
                    if (err) {
                        if (err.code === 'ER_DUP_ENTRY') {
                            console.log('Duplicate Record');
                            connection.release();
                            return res.status(400).send('Dataset title already exists');
                        } else {
                            connection.rollback(function () {
                                connection.release();
                            });
                        }
                    } else {
                        const datasetId = datasetResult.insertId;

                        for (let file of files) {
                            if (filesFromPublic) {
                                file = file.replace(/^\/?publicDatasets\//, '/');
                            }
                            connection.query('INSERT INTO file (file_loc, dataset_id) VALUES (?, ?)', [file, datasetId]);
                        }

                        // Commit transaction
                        connection.commit(function (err) {
                            if (err) {
                                connection.rollback(function () {
                                    connection.release();
                                });
                            }

                            console.log('Transaction completed successfully');
                            connection.release();
                            res.status(201).jsonp('Dataset Created.');
                        });
                    }
                });

            });
        });
    });

    if (makeItpublic) {
        try {
            let dirName = "";
            const fromPublic = false;
            if (files.length > 0) {
                dirName = path.dirname(files[0])
            }

            let userPrivateStorageDir = storageDir + username // Change this to the user's private storage path

            // Copy files from user's private storage to public dataset directory
            await copyFiles(userPrivateStorageDir, publicStorage, dirName, files, fromPublic);

        } catch (err) {
            console.error(err);
        }
    }
});

app.put('/updateDataset', async (req, res) => {

    const { title, n_cells, reference, summary, authToken, files, currentFileList } = req.body;
    const username = getUserFromToken(authToken);

    const insertList = files.filter(item => !currentFileList.includes(item));
    const deleteList = currentFileList.filter(item => !files.includes(item));

    let filesFromPublic = false;

    // Logic to Copy files from public storage to user private storage if it is a public Dataset.
    for (const file of files) {
        if (file.startsWith("publicDatasets") || file.startsWith("/publicDatasets")) {
            filesFromPublic = true;
            break;
        }
    }


    if (filesFromPublic) {
        let dirName = ""

        if (files.length > 0) {
            dirName = path.dirname(files[0])
        }

        let userPrivateStorageDir = storageDir + username // Change this to the user's private storage path

        // Copy files from user's private storage to public dataset directory
        await copyFiles("/usr/src/app/storage/", userPrivateStorageDir, dirName, files, filesFromPublic);
    }


    pool.getConnection(function (err, connection) {
        if (err) {
            console.error('Error getting DB connection:', err);
            return res.status(500).json({ message: 'Database connection error' });
        }

        connection.beginTransaction(function (err) {
            if (err) {
                console.error('Error starting transaction:', err);
                connection.release();
                return res.status(500).json({ message: 'Transaction error' });
            }

            // Run SELECT command
            connection.query('SELECT user_id FROM users WHERE username = ? LIMIT 1', [username], function (err, userRows) {
                if (err) {
                    connection.rollback(function () {
                        connection.release();
                    });
                    return res.status(500).send('Database query error');
                }

                const userId = userRows[0].user_id;

                if (!userId) {
                    res.status(400).send('User not found');
                    connection.rollback(function () {
                        connection.release();
                    });
                    return;
                }

                connection.query('SELECT dataset_id FROM dataset WHERE user_id = ? and title = ? LIMIT 1', [userId, title], function (err, datasetRows) {
                    if (err) {
                        connection.rollback(function () {
                            connection.release();
                        });
                        return res.status(500).send('Dataset query error');
                    }

                    const datasetId = datasetRows[0].dataset_id;

                    if (!datasetId) {
                        res.status(400).send('Dataset not found');
                        connection.rollback(function () {
                            connection.release();
                        });
                        return;
                    }

                    connection.query('UPDATE dataset SET n_cells=?, reference=?, summary=? WHERE dataset_id=?', [n_cells, reference, summary, datasetId], function (err, datasetResult) {
                        if (err) {
                            connection.rollback(function () {
                                connection.release();
                            });
                            return res.status(500).send('Dataset update error');
                        }
                        for (let file of insertList) {
                            if (filesFromPublic) {
                                file = file.replace(/^\/?publicDatasets\//, '/');
                            }
                            connection.query('INSERT INTO file (file_loc, dataset_id) VALUES (?, ?)', [file, datasetId]);
                        }

                        for (const file of deleteList) {
                            connection.query('DELETE FROM file WHERE file_loc=? AND dataset_id=?', [file, datasetId]);
                        }

                        // Commit transaction
                        connection.commit(function (err) {
                            if (err) {
                                connection.rollback(function () {
                                    connection.release();
                                });
                                return res.status(500).send('Dataset update error');
                            }

                            console.log('Transaction completed successfully');
                            connection.release();
                            res.status(200).jsonp('Dataset Updated.');
                        });
                    });
                });
            });
        });
    });
});

app.delete('/deleteDataset', async (req, res) => {
    const { authToken, dataset } = req.query;
    const username = getUserFromToken(authToken);

    pool.getConnection(function (err, connection) {
        if (err) {
            console.error('Error getting DB connection:', err);
            return res.status(500).send('Database connection error');
        }

        connection.beginTransaction(function (err) {
            if (err) {
                console.error('Error starting transaction:', err);
                connection.release();
                return res.status(500).send('Transaction error');
            }

            // Run SELECT command
            connection.query('SELECT user_id FROM users WHERE username = ? LIMIT 1', [username], function (err, userRows) {
                if (err) {
                    console.error('Error in SELECT query:', err);
                    connection.rollback(function () {
                        connection.release();
                    });
                    return res.status(500).send('Database query error');
                }

                const userId = userRows[0].user_id;

                if (!userId) {
                    res.status(400).send('User not found');
                    connection.rollback(function () {
                        connection.release();
                    });
                    return;
                }

                connection.query('SELECT dataset_id FROM dataset WHERE user_id = ? and title = ? LIMIT 1', [userId, dataset], function (err, datasetRows) {
                    if (err) {
                        console.error('Error in SELECT query for dataset:', err);
                        connection.rollback(function () {
                            connection.release();
                        });
                        return res.status(500).send('Dataset query error');
                    }

                    const datasetId = datasetRows[0].dataset_id;

                    if (!datasetId) {
                        res.status(400).send('Dataset not found');
                        connection.rollback(function () {
                            connection.release();
                        });
                        return;
                    }

                    connection.query('delete FROM file WHERE dataset_id=?', [datasetId], function (err, datasetResult) {
                        if (err) {
                            console.error('Error deleting files:', err);
                            connection.rollback(function () {
                                connection.release();
                            });
                            return res.status(500).send('Error deleting files');
                        }
                        connection.query('DELETE FROM dataset where dataset_id=?', [datasetId]);

                        // Commit transaction
                        connection.commit(function (err) {
                            if (err) {
                                console.error('Error committing transaction:', err);
                                connection.rollback(function () {
                                    connection.release();
                                });
                                return res.status(500).send('Transaction commit error');
                            }

                            console.log('Transaction completed successfully');
                            connection.release();
                            res.status(200).jsonp('Dataset deleted.');
                        });
                    });
                });
            });
        });
    });
});


app.post('/renameFile', async (req, res) => {
    let { oldName } = req.query;
    let { newName } = req.query;
    let { authToken } = req.query;

    oldName = oldName.replace('//', '/');
    newName = newName.replace('//', '/');

    const uname = getUserFromToken(authToken);
    if (uname == 'Unauthorized')
        return res.status(403).jsonp('Unauthorized');

    pool.query(`SELECT f.file_loc FROM aisinglecell.file f JOIN aisinglecell.dataset d ON f.dataset_id = d.dataset_id JOIN aisinglecell.users u ON d.user_id = u.user_id WHERE u.username = '${uname}' AND f.file_loc LIKE '${oldName}%';`, (err, results) => {
        if (err) {
            console.error(err);
            return res.status(500).json({ status: 500, message: 'Internal Server Error' });
        }

        if (results.length > 0) {
            return res.status(409).json({ status: 409, message: 'Directory already exists' });
        } else {
            if (oldName.includes("publicDatasets") && newName.includes("publicDatasets")) {
                fs.rename(`/usr/src/app/storage/${oldName}`, `/usr/src/app/storage/${newName}`, (err) => {
                    if (err) {
                        console.error(err);
                        return res.status(500).json({ status: 500, message: 'Internal Server Error' });
                    } else {
                        console.log('File renamed successfully!');
                        return res.status(200).jsonp('Ok');
                    }
                });
            } else {
                fs.rename(`${storageDir}${uname}/${oldName}`, `${storageDir}${uname}/${newName}`, (err) => {
                    if (err) {
                        console.error(err);
                        return res.status(500).json({ status: 500, message: 'Internal Server Error' });
                    } else {
                        console.log('File renamed successfully!');
                        return res.status(200).jsonp('Ok');
                    }
                });
            }
        }
    });
});


app.post('/download', async (req, res) => {
    const { fileList } = req.body;
    const { authToken } = req.query;
    const { pwd } = req.query;

    const username = getUserFromToken(authToken);
    if (fileList && Array.isArray(fileList)) {
        const zipName = 'files.zip';
        const output = fs.createWriteStream(zipName);
        const archive = archiver('zip');

        archive.pipe(output);

        async function appendToArchive(filePath, archivePath) {
            const fileStat = await stat(filePath);

            if (fileStat.isDirectory()) {
                // Recursively append directory contents
                const dirEntries = await fs.promises.readdir(filePath);
                for (const entry of dirEntries) {
                    const entryPath = path.join(filePath, entry);
                    const entryArchivePath = path.join(archivePath, entry);
                    await appendToArchive(entryPath, entryArchivePath);
                }
            } else {
                // Append file to archive
                archive.append(fs.createReadStream(filePath), { name: archivePath });
            }
        }

        for (const item of fileList) {
            let filePath = "";
            if (pwd && pwd.includes("publicDatasets")) {
                filePath = path.join(storageDir, item);
            } else {
                filePath = path.join(storageDir, username, item);
            }
            const archivePath = item;
            await appendToArchive(filePath, archivePath);
        }

        archive.finalize();

        output.on('close', () => {
            const zipPath = path.join(__dirname, zipName);
            const zipSize = fs.statSync(zipPath).size;
            res.setHeader('Content-disposition', 'attachment; filename=' + zipName);
            res.setHeader('Content-type', 'application/zip');
            res.setHeader('Content-length', zipSize);

            const zipstream = fs.createReadStream(zipPath);
            zipstream.pipe(res);

            // Delete the zip file after it has been sent to the client
            // fs.unlinkSync(zipPath);
        });
    } else {
        return res.status(400).jsonp('Invalid request');
    }
});

app.get('/download', async (req, res) => {
    const { fileUrl, authToken, forResultFile } = req.query;
    const { pwd } = req.query
    const username = getUserFromToken(authToken);
    let filePath = '';

    if (!fileUrl) {
        return res.status(400).jsonp('Invalid request');
    }

    if (pwd && pwd.includes("publicDatasets")) {
        filePath = path.join(storageDir, fileUrl);
    } else {
        filePath = path.join(storageDir, username, fileUrl);
    }


    const fileStat = await fs.promises.stat(filePath);

    if (fileStat.isFile()) {
        // Download file
        const filename = path.basename(fileUrl);
        const mimetype = mime.getType(filePath, { legacy: true });

        res.setHeader('Content-disposition', 'attachment; filename=' + filename);
        res.setHeader('Content-type', mimetype);

        const filestream = fs.createReadStream(filePath);
        console.log('Filename: ' + filename)
        filestream.pipe(res);
    } else if (fileStat.isDirectory()) {
        // Download folder as zip
        const folderName = path.basename(filePath);
        const archive = archiver('zip', { zlib: { level: 9 } });

        archive.directory(filePath, folderName);
        archive.pipe(res);

        res.setHeader('Content-disposition', 'attachment; filename=' + folderName + '.zip');
        res.setHeader('Content-type', 'application/zip');

        archive.finalize();
    } else {
        return res.status(400).jsonp('Invalid request');
    }
});

app.get('/fetchPreview', async (req, res) => {
    const { fileUrl, authToken, forResultFile } = req.query;
    const username = getUserFromToken(authToken);
    let filePath = '';

    if (!fileUrl) {
        return res.status(400).jsonp('Invalid request');
    }

    if (!forResultFile)
        filePath = `${storageDir}/${username}/${fileUrl}`;
    else
        filePath = `${intermediateStorage}/${fileUrl}`;

    console.log('file: ' + filePath);
    const fileStat = await fs.promises.stat(filePath);

    if (fileStat.isFile()) {
        // Read first 100 lines of the file
        const fileStream = fs.createReadStream(filePath, { encoding: 'utf8' });
        let lines = '';

        fileStream.on('data', (data) => {
            lines += data;

            // Check if 100 lines have been read
            if (lines.split('\n').length >= 20) {
                fileStream.destroy();
            }
        });

        fileStream.on('close', () => {
            res.status(200).send(lines);
        });

        fileStream.on('error', (error) => {
            console.log('Error reading file: ' + error);
            res.status(500).jsonp('Error reading file');
        });
    } else {
        return res.status(400).jsonp('Invalid request');
    }
});

app.delete('/deleteFiles', async (req, res) => {
    const { fileList } = req.body;
    const { authToken } = req.query;
    const { pwd } = req.query
    const uname = getUserFromToken(authToken);
    if (uname === 'Unauthorized') {
        return res.status(403).json('Unauthorized');
    }
    if (fileList && Array.isArray(fileList)) {
        let successCount = 0;
        let errorCount = 0;
        let failFlag = false;
        for (const file of fileList) {
            try {
                const [rows, fields] = await pool.promise().execute(`SELECT f.file_loc FROM aisinglecell.file f JOIN aisinglecell.dataset d ON f.dataset_id = d.dataset_id JOIN aisinglecell.users u ON d.user_id = u.user_id WHERE u.username = '${uname}' AND f.file_loc like '${file.replace("'", "''")}%';`);
                console.log('Count: ' + rows.length);
                if (rows.length > 0) {
                    res.status(401).json({ message: 'File(s) being used by datasets.' });
                    return;
                } else {
                    let filePath = ""
                    if (pwd.includes("publicDatasets")) {
                        filePath = `${storageDir}${file}`;
                    } else {
                        filePath = `${storageDir}${uname}/${file}`;
                    }
                    console.log(filePath);
                    try {
                        if (fs.existsSync(filePath)) {
                            fs.removeSync(filePath);
                            // file or folder removed
                            successCount++;
                        } else {
                            console.log('Error Deleting: ' + filePath);
                            errorCount++;
                        }
                    } catch (err) {
                        console.error(err);
                        errorCount++;
                    }
                }
            } catch (err) {
                console.error(err);
                res.status(500).json({ message: 'Internal Server Error' });
                return;
            }
        }
        return res.jsonp({ success: successCount, errorCount: errorCount });
    } else {
        return res.status(400).jsonp('Invalid request');
    }
});



const formatDate = (dateString) => {
    const date = new Date(dateString);
    const options = {
        day: 'numeric',
        month: 'short',
        year: '2-digit',
        hour: 'numeric',
        minute: 'numeric',
        hour12: false
    };
    const formattedDate = date.toLocaleDateString('en-GB', options);
    return formattedDate;
}

app.get('/getDirContents', async (req, res) => {
    try {
        console.log(`HOSTURL: ${process.env.HOST_URL}`);
        const { dirPath, authToken, usingFor } = req.query;

        let uid = getUserFromToken(authToken);
        if (uid == "Unauthorized") {
            return res.status(403).jsonp(uid);
        }

        subdir = req.query.subdir;
        var directoryPath = ""


        var directoryPath = path.join(storageDir + uid + "/" + dirPath + "/");

        if (subdir != undefined)
            directoryPath = path.join(storageDir + uid + "/", subdir);

        if (dirPath == "publicDatasets") {
            directoryPath = publicStorage;
        }

        if (dirPath.includes("publicDatasets/")) {
            directoryPath = "/usr/src/app/storage/" + dirPath;
        }

        // Check if the directory exists
        if (!fs.existsSync(directoryPath)) {
            // Create the directory if it doesn't exist
            fs.mkdirSync(directoryPath, { recursive: true });
            console.log(`Directory "${directoryPath}" created successfully.`);
        } else {
            console.log(`Directory "${directoryPath}" already exists.`);
        }

        const directoryContents = fs.readdirSync(directoryPath);
        const dirList = [];
        const fileList = [];


        directoryContents.forEach((item) => {
            const itemPath = `${directoryPath}/${item}`;
            const itemStats = fs.statSync(itemPath);


            if (itemStats.isDirectory() == true)
                dirList.push({ "name": item, "created": formatDate(itemStats.birthtime) });
            else {
                let dotIndex = item.lastIndexOf('.');
                fileList.push({ "name": item, "created": formatDate(itemStats.birthtime), "type": (dotIndex != -1 ? item.substring(dotIndex + 1).toUpperCase() + " " : "") });
            }
        });

        return res.json({ 'Directories': dirList, 'Files': fileList });

    }
    // uid = req.session.username;
    catch (e) {
        console.log('Errordsd: ' + e);
        // return res.status(400).jsonp('Invalid request');
    }

});

app.post('/upload', async (req, res) => {
    let { uploadDir, authToken, publicDatasetFlag } = req.query;
    let username = getUserFromToken(authToken);

    let destDir = publicDatasetFlag === "true" ? "./storage/" + uploadDir : "./storage/" + username + uploadDir;

    let tempDir = './uploads';

    let storage = multer.diskStorage({
        destination: (req, file, cb) => {
            cb(null, tempDir);
        },
        filename: (req, file, cb) => {
            console.log(file.originalname);
            cb(null, file.originalname);
        },
    });

    let uploadFiles = multer({
        storage: storage,
        // limits: { fileSize: FILE_UPLOAD_MAX_SIZE },
    }).array("files", 10);

    let uploadFunction = util.promisify(uploadFiles);

    try {
        await uploadFunction(req, res);

        // Move uploaded files to storage directory
        let files = req.files;
        for (let i = 0; i < files.length; i++) {
            let file = files[i];
            let tempFilePath = path.join(tempDir, file.filename);
            let destFilePath = path.join(destDir, file.originalname);

            console.log(`Tempstorage: ${tempFilePath}, DestinationL ${destFilePath}`);
            if (!fs.existsSync(path.join(destDir))) {
                fs.mkdirSync(path.join(destDir), { recursive: true });
            }

            fs.copyFileSync(tempFilePath, destFilePath);
            fs.unlinkSync(tempFilePath);
        }

        res.status(200).json({ message: 'File uploaded successfully' });
    } catch (error) {
        res.status(500).json({ message: 'Failed to upload file', error });
    }
});


app.post('/createNewFolder', (req, res) => {
    const { pwd, folderName, authToken } = req.query;
    const username = getUserFromToken(authToken);
    let folderPath = ""
    if (pwd.includes("publicDatasets")) {
        folderPath = `/usr/src/app/storage/${pwd}/${folderName}`;
    } else {
        folderPath = `${storageDir}/${username}/${pwd}/${folderName}`;
    }
    if (fs.existsSync(folderPath)) {
        res.status(400).jsonp('Folder already exists');
        return;
    }
    try {
        fs.promises.mkdir(folderPath);
        res.status(201).jsonp('Folder created')
    }
    catch (err) {
        res.status(404).jsonp('Bad root folder: ' + err);
    }
})

const { exec } = require('child_process');

app.get('/getStorageDetails', async (req, res) => {
    const sizeRegex = /^(\d+(\.\d+)?)\s*([KMG]B?)$/i;
    try {
        const { authToken } = req.query;

        let username = getUserFromToken(authToken);
        if (username === "Unauthorized") {
            return res.status(403).jsonp(uid);
        }

        const cmd = `du -sh ${storageDir}/${username}`;
        exec(cmd, (err, stdout, stderr) => {
            if (err) {
                console.error(err);
                return res.status(500).json({ error: 'Internal server error' });
            }
            if (stderr) {
                console.error(stderr);
                return res.status(500).json({ error: 'Internal server error' });
            }
            const [size, folder] = stdout.split('\t');
            const match = size.match(sizeRegex);
            const [_, value, __, unit] = match;
            const gigabytes = parseFloat(value) / ({ K: 1024 * 1024, M: 1024, G: 1 }[unit.toUpperCase()]);

            console.log(`disk utilization: ${gigabytes} GB, folder: ${folder}`);

            return res.json({ used: gigabytes, allowed: storageAllowance });
        });
    } catch (e) {
        console.log('Error in getting Storage usage: ' + e);
        return res.status(400).jsonp('Invalid request');
    }
});

// Route to get datasets and files for a specific user
app.get('/preview/datasets', (req, res) => {

    const { authToken } = req.query;

    const username = getUserFromToken(authToken);

    if (username == "Unauthorized") {
        return res.status(403).jsonp(username);
    }

    // Get user ID based on username
    const userQuery = `SELECT user_id FROM users WHERE username = '${username}'`;

    pool.query(userQuery, (err, userResult) => {
        if (err) {
            console.error('Database query error:', err);
            return res.status(500).json({ message: 'Internal Server Error' });
        }

        if (userResult.length === 0) {
            res.status(404).send(`User '${username}' not found`);
        } else {
            const userID = userResult[0].user_id;

            // Get datasets and files for the specified user
            const datasetsQuery = `
          SELECT dataset.dataset_id, dataset.title, dataset.n_cells, dataset.reference, dataset.summary, file.file_id, file.file_loc, SUBSTRING_INDEX(SUBSTRING_INDEX(file.file_loc, '/', 2), '/', -1) AS direc
          FROM dataset
          JOIN file ON dataset.dataset_id = file.dataset_id
          WHERE dataset.user_id = ${userID}
        `;

            pool.query(datasetsQuery, (err, datasetsResult) => {
                if (err) {
                    console.error('Database query error:', err);
                    return res.status(500).json({ message: 'Internal Server Error' });
                }

                const datasets = {};

                datasetsResult.forEach(row => {
                    const { dataset_id, title, n_cells, reference, summary, file_id, file_loc, direc } = row;
                    if (!datasets[dataset_id]) {
                        datasets[dataset_id] = {
                            title,
                            n_cells,
                            reference,
                            summary,
                            files: [],
                            direc,
                            dataset_id
                        };
                    }
                    datasets[dataset_id].files.push({
                        file_id,
                        file_loc
                    });
                });

                res.json(datasets);
            });
        }
    });
});


// Define API endpoint
app.get('/api/tools/leftnav', function (req, res) {
    // Query the category and filter tables and group the filters by category
    const sql = 'SELECT c.id AS category_id, c.name AS category_name, ' +
        'JSON_ARRAYAGG(f.name) AS filters ' +
        'FROM categories c ' +
        'LEFT JOIN filters f ON c.id = f.category_id ' +
        'GROUP BY c.id ' +
        'ORDER BY c.id ASC';
    pool.query(sql, function (error, results, fields) {
        if (error) {
            console.log(error);
            res.status(500).send('Internal server error');
        } else {
            res.json(results);
        }
    });
});

app.post('/createTask', (req, res) => {
    const { taskTitle, taskId, method, authToken, outputPath } = req.body;
    const username = getUserFromToken(authToken);

    pool.getConnection(function (err, connection) {
        if (err) {
            console.error('Error getting DB connection:', err);
            return res.status(500).json({ message: 'Database connection error' });
        }

        connection.beginTransaction(function (err) {
            if (err) {
                console.error('Error starting transaction:', err);
                connection.release();
                return res.status(500).json({ message: 'Transaction error' });
            }

            connection.query('SELECT user_id FROM users WHERE username = ? LIMIT 1', [username], function (err, userRows) {
                if (err) {
                    console.error('Error in SELECT query:', err);
                    connection.rollback(function () {
                        connection.release();
                    });
                    return res.status(500).json({ message: 'Database query error' });
                }

                const userId = userRows[0].user_id;

                if (!userId) {
                    res.status(400).send('User not found');
                    connection.rollback(function () {
                        connection.release();
                    });
                    return;
                }

                const date = new Date();
                const timestamp = Date.UTC(date.getUTCFullYear(), date.getUTCMonth(), date.getUTCDate(), date.getUTCHours(), date.getUTCMinutes(), date.getUTCSeconds(), date.getUTCMilliseconds());
                connection.query('INSERT INTO task (task_title, task_id, user_id, tool, results_path, created_datetime) VALUES (?,?, ?, ?, ?, ?)', [taskTitle, taskId, userId, method, outputPath, timestamp], function (err, taskResult) {
                    if (err) {
                        console.error('Error in INSERT query:', err);
                        connection.rollback(function () {
                            connection.release();
                        });
                        return res.status(500).json({ message: 'Database insertion error' });
                    } else {
                        // Commit transaction
                        connection.commit(function (err) {
                            if (err) {
                                console.error('Error committing transaction:', err);
                                connection.rollback(function () {
                                    connection.release();
                                });
                                return res.status(500).json({ message: 'Transaction commit error' });
                            }

                            console.log('Transaction completed successfully');
                            connection.release();
                            res.status(201).jsonp('Task Created.');
                        });
                    }
                });
            });
        });
    });
});

app.put('/updateTaskStatus', (req, res) => {
    const { taskIds, status } = req.body;
    const taskIdsArr = taskIds.split(',');

    pool.getConnection(function (err, connection) {
        if (err) {
            console.error('Error getting DB connection:', err);
            return res.status(500).json({ message: 'Database connection error' });
        }

        connection.beginTransaction(function (err) {
            if (err) {
                console.error('Error starting transaction:', err);
                connection.release();
                return res.status(500).json({ message: 'Transaction error' });
            }

            const date = new Date();
            const timestamp = Date.UTC(date.getUTCFullYear(), date.getUTCMonth(), date.getUTCDate(), date.getUTCHours(), date.getUTCMinutes(), date.getUTCSeconds(), date.getUTCMilliseconds());
            connection.query('UPDATE task SET status = ?, finish_datetime = ? WHERE task_id IN (?)', [status, timestamp, taskIdsArr], function (err, taskResult) {
                if (err) {
                    console.error('Error in UPDATE query:', err);
                    connection.rollback(function () {
                        connection.release();
                    });
                    return res.status(500).json({ message: 'Database update error' });
                } else {
                    // Commit transaction
                    connection.commit(function (err) {
                        if (err) {
                            console.error('Error committing transaction:', err);
                            connection.rollback(function () {
                                connection.release();
                            });
                            return res.status(500).json({ message: 'Transaction commit error' });
                        }

                        console.log('Transaction completed successfully');
                        connection.release();
                        res.status(200).jsonp('Task status updated.');
                    });
                }
            });
        });
    });
});

app.get('/getTasks', (req, res) => {
    const { authToken } = req.query;
    const username = getUserFromToken(authToken);

    pool.getConnection(function (err, connection) {
        if (err) {
            console.error('Error getting DB connection:', err);
            return res.status(500).json({ message: 'Database connection error' });
        }

        connection.beginTransaction(function (err) {
            if (err) {
                console.error('Error starting transaction:', err);
                connection.release();
                return res.status(500).json({ message: 'Transaction error' });
            }
            connection.query('SELECT user_id FROM users WHERE username = ? LIMIT 1', [username], function (err, userRows) {
                if (err) {
                    console.error('Database query error:', err);
                    return res.status(500).json({ message: 'Internal Server Error' });
                }

                const userId = userRows[0].user_id;

                if (!userId) {
                    res.status(400).send('User not found');
                    connection.rollback(function () {
                        connection.release();
                    });
                    return;
                }

                connection.query('SELECT task_title, task_id, results_path, tool, status, created_datetime, finish_datetime FROM task WHERE user_id = ?', [userId], function (err, rows) {
                    if (err) {
                        console.error('Error committing transaction:', err);
                        connection.rollback(function () {
                            connection.release();
                        });
                        return res.status(500).json({ message: 'Transaction commit error' });
                    } else {
                        connection.release();
                        res.json(rows);
                    }
                });
            });
        });
    });
});

//New endpoint to check if the hash for a file exists or not
app.get('/mongoDB/api/file-exists', async (req, res) => {
    const hash = req.query.hash;
    console.log(hash)

    const client = new MongoClient(mongoUrl, { useUnifiedTopology: true });

    // Connect to the MongoDB server
    await client.connect();

    const db = client.db(dbName);
    const collection = db.collection(userDatasetsCollection);

    // Query MongoDB to check if the hash exists
    const result = await collection.count({
        fileHashes: { $in: [hash] }
    })


    if (result >= 1) {
        res.status(200).json({ exists: true });
    } else {
        res.status(200).json({ exists: false });
    }

    client.close()
});


// Connect to MongoDB and retrieve options
app.get('/mongoDB/api/options', async (req, res) => {
    try {
        const client = new MongoClient(mongoUrl, { useUnifiedTopology: true });

        // Connect to the MongoDB server
        await client.connect();

        const db = client.db(dbName);
        const collection = db.collection(optionsCollectionName);

        // Define the unique compound index on 'field' and 'name'
        await collection.createIndex({ field: 1, name: 1 }, { unique: true });

        // Use the aggregation framework to group options by field
        const pipeline = [
            {
                $group: {
                    _id: '$field',
                    options: { $addToSet: { name: '$name', abbreviation: '$abbreviation' } },
                },
            },
        ];

        const result = await collection.aggregate(pipeline).toArray();

        // Transform the result into an object with field names as keys
        const optionsByField = {};
        result.forEach((item) => {
            optionsByField[item._id] = item.options;
        });

        // Close the MongoDB connection
        client.close();

        // Return the options as a JSON response
        res.status(200).json(optionsByField);
    } catch (err) {
        console.error('Error:', err);
        res.status(500).json({ error: 'Internal Server Error' });
    }
});

app.post('/mongoDB/api/submitDatasetMetadata', async (req, res) => {
    const client = new MongoClient(mongoUrl);

    try {
        const formData = req.body; // This assumes you have middleware to parse JSON in the request body
        const makeItpublic = formData.makeItpublic;
        let files = formData.files;
        let authToken = formData.userId;

        if (authToken) {
            let username = getUserFromToken(authToken);
            formData.userId = username;
        }


        // Connect to the MongoDB server
        await client.connect();
        const db = client.db(dbName);

        const collection = formData.flow == 'upload' ? db.collection(datasetCollection) : db.collection(userDatasetsCollection);



        // Check if a document with the provided Id already exists
        const existingDocument = await collection.findOne({ Id: formData.Id });

        if (existingDocument) {
            console.log('Document with Id already exists:', formData.Id);
            res.status(400).json({ error: 'Document with the provided Id already exists' });
        } else {
            // Document with the provided Id does not exist, proceed with insertion
            await collection.insertOne(formData);
            console.log('Form data submitted successfully');

            if (makeItpublic) {
                console.log("Transfering files from local to public folder");
                try {
                    let dirName = "";
                    const fromPublic = false;
                    if (files && files.length > 0) {
                        dirName = path.dirname(files[0])
                    }

                    let userPrivateStorageDir = storageDir + username // Change this to the user's private storage path

                    // Copy files from user's private storage to public dataset directory
                    await copyFiles(userPrivateStorageDir, publicStorage, dirName, files, fromPublic);

                } catch (err) {
                    console.error(err);
                }
            }
            res.status(200).json({ message: 'Form data submitted successfully' });
        }
    } catch (err) {
        console.error('Error:', err);
        res.status(500).json({ error: 'Internal Server Error' });
    } finally {
        // Ensure the client will close when you finish/error
        await client.close();
    }
});


app.post('/mongoDB/api/submitTaskMetadata', async (req, res) => {
    const client = new MongoClient(mongoUrl);

    try {
        await client.connect();
        const db = client.db(dbName);
        const collection = db.collection(taskBuilderCollectionName);

        let documents = req.body;
        // Ensure documents is always an array for consistency
        if (!Array.isArray(documents)) {
            documents = [documents];
        }

        const insertionResults = [];
        for (const formData of documents) {
            // Check if a document with the provided Id already exists
            const existingDocument = await collection.findOne({ Id: formData.Id });

            if (existingDocument) {
                console.log('Document with Id already exists:', formData.Id);
                insertionResults.push({
                    Id: formData.Id,
                    status: 'error',
                    message: 'Document with the provided Id already exists',
                });
            } else {
                // Document with the provided Id does not exist, proceed with insertion
                await collection.insertOne(formData);
                console.log('Form data submitted successfully for Id:', formData.Id);
                insertionResults.push({
                    Id: formData.Id,
                    status: 'success',
                    message: 'Form data submitted successfully',
                });
            }
        }

        // If handling multiple documents, you might want to aggregate results and respond accordingly
        if (insertionResults.length > 1) {
            // Respond with the aggregated results for multiple documents
            res.status(200).json(insertionResults);
        } else if (insertionResults.length === 1) {
            // For a single document, you can respond with the single result
            const result = insertionResults[0];
            if (result.status === 'success') {
                res.status(200).json({ message: result.message });
            } else {
                res.status(400).json({ error: result.message });
            }
        } else {
            // No documents were processed
            res.status(400).json({ error: 'No documents were submitted' });
        }
    } catch (err) {
        console.error('Error:', err);
        res.status(500).json({ error: 'Internal Server Error' });
    } finally {
        await client.close();
    }
});




// API endpoint to get datasets
app.get('/mongoDB/api/getDatasets', async (req, res) => {
    const client = new MongoClient(mongoUrl);

    try {
        // Connect to the MongoDB server
        await client.connect();
        const db = client.db(dbName);
        const collection = db.collection(datasetCollection);

        // Fetching all documents from the collection
        const datasets = await collection.find({}).toArray();

        // Sending the datasets as a JSON response
        res.status(200).json(datasets);
    } catch (err) {
        console.error('Error:', err);
        res.status(500).json({ error: 'Internal Server Error' });
    } finally {
        // Ensure the client will close when you finish/error
        await client.close();
    }
});
// Define a route to handle adding a new option to MongoDB
app.post('/mongoDB/api/addNewOption', async (req, res) => {
    const { field, name, username } = req.body;

    // Create the document object with the specified format
    const newOption = {
        field: field,
        name: name,
        username: username
    };
    try {
        const client = new MongoClient(mongoUrl, { useUnifiedTopology: true });
        await client.connect();

        const db = client.db(dbName);
        const collection = db.collection(optionsCollectionName);

        // Define the unique compound index on 'field' and 'name'
        await collection.createIndex({ field: 1, name: 1 }, { unique: true });

        // Insert the new option into the collection
        const insertResult = await collection.insertOne(newOption);

        client.close();

        res.status(200).json({
            message: `New option "${name}" added to MongoDB for field "${field}"`,
            insertedId: insertResult.insertedId,
        });
    } catch (error) {
        console.error('Error adding new option to MongoDB:', error);
        res.status(500).json({ error: 'Internal Server Error' });
    }
});

// Connect to MongoDB and retrieve options
app.get('/mongoDB/api/groupedUserOptions', async (req, res) => {
    try {
        const client = new MongoClient(mongoUrl, { useUnifiedTopology: true });

        // Connect to the MongoDB server
        await client.connect();

        const db = client.db(dbName);
        const collection = db.collection(optionsCollectionName);

        const username = req.query.username;
        const isAdmin = req.query.isAdmin;

        // Define the match stage of the aggregation pipeline
        const matchStage = isAdmin === 'true' ? {} : { username: username };

        // Aggregation pipeline stages
        const pipeline = [
            { $match: matchStage },
            {
                $group: {
                    _id: '$field',
                    options: { $addToSet: { _id: '$_id', name: '$name', username: '$username', abbreviation: '$abbreviation' } },
                },
            },
        ];

        const result = await collection.aggregate(pipeline).toArray();

        // Transform the result into an object with field names as keys
        const optionsByField = {};
        result.forEach((item) => {
            optionsByField[item._id] = item.options;
        });

        // Close the MongoDB connection
        client.close();

        // Return the options as a JSON response
        res.status(200).json(optionsByField);
    } catch (err) {
        console.error('Error:', err);
        res.status(500).json({ error: 'Internal Server Error' });
    }
});

// Define a DELETE route to delete selected options
app.delete('/mongoDB/api/deleteOptions', async (req, res) => {
    try {
        const optionIds = req.body.optionIds; // Assuming the request body contains an array of option IDs

        const client = new MongoClient(mongoUrl, { useUnifiedTopology: true });

        // Connect to the MongoDB server
        await client.connect();

        const db = client.db(dbName);
        const collection = db.collection(optionsCollectionName);

        // Convert optionIds to MongoDB ObjectIDs
        const objectIds = optionIds.map(id => new ObjectId(id));

        // Delete the options with the specified ObjectIDs
        const deleteResult = await collection.deleteMany({ _id: { $in: objectIds } });

        client.close();

        if (deleteResult.deletedCount > 0) {
            res.status(200).json({ message: 'Options deleted successfully' });
        } else {
            res.status(404).json({ message: 'Options not found' });
        }
    } catch (error) {
        console.error('Error deleting options:', error);
        res.status(500).json({ message: 'Internal server error' });
    }
});


// Define a route to handle adding a new option for Task field to MongoDB
app.post('/mongoDB/api/addTaskOption', async (req, res) => {
    const { field, name, username, abbreviation } = req.body;

    // Create the document object with the specified format
    const newOption = {
        field: field,
        name: name,
        username: username,
        abbreviation: abbreviation
    };
    try {
        const client = new MongoClient(mongoUrl, { useUnifiedTopology: true });
        await client.connect();

        const db = client.db(dbName);
        const collection = db.collection(optionsCollectionName);

        // Define the unique compound index on 'field' and 'name'
        await collection.createIndex({ field: 1, name: 1 }, { unique: true });

        // Insert the new option into the collection
        const insertResult = await collection.insertOne(newOption);

        client.close();

        res.status(200).json({
            message: `New option "${name}" added to MongoDB for field "${field}"`,
            insertedId: insertResult.insertedId,
        });
    } catch (error) {
        console.error('Error adding new option to MongoDB:', error);
        res.status(500).json({ error: 'Internal Server Error' });
    }
});


//API to move files from one folder to another
app.post('/api/move-files', (req, res) => {
    const { newDirectoryPath, jwtToken } = req.body;
    const username = getUserFromToken(jwtToken);
    let destinationPath = ""
    if (username) {
        destinationPath = `${storageDir}/${username}/${newDirectoryPath}`;
    }
    let sourcePath = `${storageDir}/tempStorage`;

    if (!fs.existsSync(destinationPath)) {
        fs.mkdirSync(destinationPath, { recursive: true });
    }

    const files = fs.readdirSync(sourcePath);

    files.forEach((filename) => {
        const sourcePathFile = path.join(sourcePath, filename);
        const destinationPathFile = path.join(destinationPath, filename);

        fs.renameSync(sourcePathFile, destinationPathFile);
    });

    res.sendStatus(200);
});


app.delete('/api/storage/delete-file', (req, res) => {

    try {
        const { fileName, authToken, newDirectoryPath } = req.query;

        const uname = getUserFromToken(authToken);
        if (uname == 'Unauthorized')
            return res.status(403).jsonp('Unauthorized');

        let filepath = `${storageDir}/${newDirectoryPath}/${fileName}`;


        fs.unlink(filepath, (err) => {
            if (err) {
                console.error("Error deleting file:", err);
                return res.status(500).send('Error deleting file');
            }
            res.send('File deleted successfully');
        });
    } catch (error) {
        console.error('Error deleting a file:', error);
        res.status(500).json({ error: 'Internal Server Error' });
    }
});

app.post('/api/storage/renameFile', async (req, res) => {

    try {

        let { oldName } = req.query;
        let { newName } = req.query;
        let { authToken } = req.query;

        const uname = getUserFromToken(authToken);

        if (uname == 'Unauthorized')
            return res.status(403).jsonp('Unauthorized');

        fs.rename(`${storageDir}${oldName}`, `${storageDir}${newName}`, (err) => {
            if (err) {
                console.error(err);
                res.status(500).json({ error: 'Internal Server Error' });
            } else {
                console.log('File renamed successfully!');
                return res.status(200).jsonp('File renamed successfully!');
            }
        });

    } catch (error) {
        console.error('Error renaming a file:', error);
        res.status(500).json({ error: 'Internal Server Error' });
    }
});

// Fetch facets and paginated results
app.post('/api/datasets/search', async (req, res) => {
    let client;
    try {
        client = new MongoClient(mongoUrl);
        await client.connect();

        const db = client.db(dbName);
        const datasetType = req.query.datasetType; // New parameter
        const collection = datasetType === "myDatasets" ? db.collection(userDatasetsCollection) : db.collection(datasetCollection);


        const page = parseInt(req.query.page, 10) || 1;
        const pageSize = parseInt(req.query.pageSize, 10) || 10;
        let globalSearchQuery = req.query.q;
        const filters = req.body.filters;

        //Update this field accordingly whenever you add a new facet 
        const fieldsWithLabel = ['Species', 'Anatomical Entity', 'Organ Part', 'Selected Cell Types', 'Disease Status (Specimen)', 'Disease Status (Donor)'];

        let matchConditions = [];

        // Add the global search query to the match conditions
        if (globalSearchQuery) {
            matchConditions.push({
                $or: [
                    { 'Species.label': { $regex: globalSearchQuery, $options: 'i' } },
                    { 'Author': { $regex: globalSearchQuery, $options: 'i' } },
                    { 'Anatomical Entity.label': { $regex: globalSearchQuery, $options: 'i' } },
                    { 'Organ Part.label': { $regex: globalSearchQuery, $options: 'i' } },
                    { 'Selected Cell Types.label': { $regex: globalSearchQuery, $options: 'i' } },
                    { 'Disease Status (Specimen).label': { $regex: globalSearchQuery, $options: 'i' } },
                    { 'Disease Status (Donor).label': { $regex: globalSearchQuery, $options: 'i' } },
                ],
            });
        }

        // Apply additional filters
        if (filters) {
            Object.keys(filters).forEach((filterCategory) => {
                const filterValue = filters[filterCategory];
                if (Array.isArray(filterValue) && filterValue.length > 0) {
                    let condition = {};

                    // Check if the filter category should use the 'label' property
                    if (fieldsWithLabel.includes(filterCategory)) {
                        condition[`${filterCategory}.label`] = { $in: filterValue };
                    } else {
                        // Directly use the filter category for other fields
                        condition[filterCategory] = { $in: filterValue };
                    }

                    // Add this condition to the matchConditions array
                    matchConditions.push(condition);
                }
            });
        }

        // Construct the final match stage using $and, only if there are multiple conditions
        let matchStage = {};
        if (matchConditions.length > 0) {
            matchStage = matchConditions.length > 1 ? { $and: matchConditions } : matchConditions[0];
        }

        // Define the pipeline for facets
        const facetsPipeline = [
            { $match: matchStage },
            {
                $facet: {
                    'Species': [
                        { $group: { _id: '$Species.label', count: { $sum: 1 } } }, { $sort: { count: -1 } }
                    ],
                    'Author': [
                        { $group: { _id: '$Author', count: { $sum: 1 } } }, { $sort: { count: -1 } }
                    ],
                    'Anatomical Entity': [
                        { $group: { _id: '$Anatomical Entity.label', count: { $sum: 1 } } }, { $sort: { count: -1 } }
                    ],
                    'Organ Part': [
                        { $group: { _id: '$Organ Part.label', count: { $sum: 1 } } }, { $sort: { count: -1 } }
                    ],
                    'Selected Cell Types': [
                        { $group: { _id: '$Selected Cell Types.label', count: { $sum: 1 } } }, { $sort: { count: -1 } }
                    ],
                    'Disease Status (Specimen)': [
                        { $group: { _id: '$Disease Status (Specimen).label', count: { $sum: 1 } } }, { $sort: { count: -1 } }
                    ],
                    'Disease Status (Donor)': [
                        { $group: { _id: '$Disease Status (Donor).label', count: { $sum: 1 } } }, { $sort: { count: -1 } }
                    ],
                    // ... add other facets here
                }
            },
        ];

        // Get the facets
        const facetsResult = await collection.aggregate(facetsPipeline).toArray();

        // Pagination: Get total count for the query
        const totalCount = await collection.countDocuments(matchStage);

        // Build the pipeline for search results with pagination
        const searchResultsPipeline = [
            { $match: matchStage },
            // { $project: { Cells: 0, Genes: 0, QC_Plots: 0, cell_metadata_obs:0, gene_metadata:0, layers:0, inputFiles:0, adata_path:0 } }, // Excluding fields
            { $skip: (page - 1) * pageSize },
            { $limit: pageSize },
        ];

        // Get the paginated search results
        const searchResults = await collection.aggregate(searchResultsPipeline).toArray();

        res.json({
            facets: facetsResult[0],
            results: searchResults,
            pagination: {
                page,
                pageSize,
                pageCount: Math.ceil(totalCount / pageSize),
                totalCount
            }
        });
    } catch (error) {
        console.error('Search failed:', error);
        res.status(500).send('An error occurred while searching.');
    } finally {
        // Ensure the MongoDB client is always closed, even if an error occurs
        if (client) {
            await client.close();
        }
    }
});


app.post('/api/tasks/search', async (req, res) => {
    let client;
    try {
        client = new MongoClient(mongoUrl, { useNewUrlParser: true, useUnifiedTopology: true });
        await client.connect();

        const db = client.db(dbName);
        const tasksCollection = db.collection(taskBuilderCollectionName);

        const page = parseInt(req.query.page, 10) || 1;
        const taskType = req.query.task_type; // Extract task_type from query parameters
        const pageSize = parseInt(req.query.pageSize, 10) || 10;
        const globalSearchQuery = req.query.q || '';
        const filters = req.body.filters;

        //Update this field accordingly whenever you add a new facet 
        const fieldsWithLabel = ['Species', 'Anatomical Entity', 'Organ Part', 'Selected Cell Types', 'Disease Status (Specimen)', 'Disease Status (Donor)'];

        if (!taskType) {
            return res.status(400).send('task_type is required');
        }
        // Initialize matchConditions to include the taskType filter
        let matchConditions = [{ 'TaskType.label': taskType }];

        // Include global search query in match conditions
        if (globalSearchQuery) {
            matchConditions.push({
                $or: [
                    { 'TaskType.label': { $regex: globalSearchQuery, $options: 'i' } },
                    { 'datasetDetails.Species.label': { $regex: globalSearchQuery, $options: 'i' } },
                    { 'datasetDetails.Author': { $regex: globalSearchQuery, $options: 'i' } },
                    { 'datasetDetails.Anatomical Entity.label': { $regex: globalSearchQuery, $options: 'i' } },
                    { 'datasetDetails.Organ Part.label': { $regex: globalSearchQuery, $options: 'i' } },
                    { 'datasetDetails.Selected Cell Types.label': { $regex: globalSearchQuery, $options: 'i' } },
                    { 'datasetDetails.Disease Status (Specimen).label': { $regex: globalSearchQuery, $options: 'i' } },
                    { 'datasetDetails.Disease Status (Donor).label': { $regex: globalSearchQuery, $options: 'i' } },
                ]
            });
        }


        // Apply additional filters
        if (filters && Object.keys(filters).length > 0) {
            Object.keys(filters).forEach((filterCategory) => {
                const filterValue = filters[filterCategory];
                if (Array.isArray(filterValue) && filterValue.length > 0) {
                    let condition = {};
                    if (fieldsWithLabel.includes(filterCategory)) {
                        condition[`datasetDetails.${filterCategory}.label`] = { $in: filterValue };
                    } else {
                        condition[`datasetDetails.${filterCategory}`] = { $in: filterValue };
                    }
                    matchConditions.push(condition);
                }
            });
        }

        // Construct the final match stage
        let matchStage = {};
        if (matchConditions.length > 1) {
            matchStage = { $and: matchConditions };
        } else if (matchConditions.length === 1) {
            matchStage = matchConditions[0];
        }

        const basePipeline = [
            {
                $lookup: {
                    from: datasetCollection,
                    localField: "DatasetId",
                    foreignField: "Id",
                    as: "datasetDetails"
                }
            },
            { $unwind: "$datasetDetails" },
            {
                $match: matchStage
            }
        ];

        // Define the pipeline for counting the total number of matching tasks
        const countPipeline = [
            ...basePipeline,
            { $count: "total" }
        ];

        // Execute the count pipeline
        const countResults = await tasksCollection.aggregate(countPipeline).toArray();
        const totalTasksCount = countResults.length > 0 ? countResults[0].total : 0;


        const facetAndDocumentsPipeline = [
            {
                $facet: {
                    // Each facet is a direct property of the `$facet` object
                    Species: [
                        { $group: { _id: '$datasetDetails.Species.label', count: { $sum: 1 } } },
                        { $sort: { count: -1 } }
                    ],
                    Author: [
                        { $group: { _id: '$datasetDetails.Author', count: { $sum: 1 } } },
                        { $sort: { count: -1 } }
                    ],
                    'Anatomical Entity': [
                        { $group: { _id: '$datasetDetails.Anatomical Entity.label', count: { $sum: 1 } } },
                        { $sort: { count: -1 } }
                    ],
                    'Organ Part': [
                        { $group: { _id: '$datasetDetails.Organ Part.label', count: { $sum: 1 } } },
                        { $sort: { count: -1 } }
                    ],
                    'Selected Cell Types': [
                        { $group: { _id: '$datasetDetails.Selected Cell Types.label', count: { $sum: 1 } } },
                        { $sort: { count: -1 } }
                    ],
                    'Disease Status (Specimen)': [
                        { $group: { _id: '$datasetDetails.Disease Status (Specimen).label', count: { $sum: 1 } } },
                        { $sort: { count: -1 } }
                    ],
                    'Disease Status (Donor)': [
                        { $group: { _id: '$datasetDetails.Disease Status (Donor).label', count: { $sum: 1 } } },
                        { $sort: { count: -1 } }
                    ],
                    // Documents facet for pagination
                    documents: [
                        { $skip: (page - 1) * pageSize },
                        { $limit: pageSize },
                        {
                            $project: {
                                Title: "$datasetDetails.Title",
                                TaskId: "$Id",
                                TaskType: "$TaskType.label",
                                Species: "$datasetDetails.Species.label",
                                'Organ Part': "$datasetDetails.Organ Part.label",
                                'Cell Count Estimate': "$datasetDetails.Cell Count Estimate",
                                'Development Stage': "$datasetDetails.Development Stage",
                                'Anatomical Entity': "$datasetDetails.Anatomical Entity.label",
                                'Disease Status (Donor)': "$datasetDetails.Disease Status (Donor).label",
                                Author: "$datasetDetails.Author",
                                TaskLabel: "$TaskLabel.label",
                                'Source': "$datasetDetails.Source",
                                'Submission Date': "$datasetDetails.Submission Date",
                            }
                        }
                    ]
                }
            }
        ];


        const finalPipeline = basePipeline.concat(facetAndDocumentsPipeline);

        const aggregatedResults = await tasksCollection.aggregate(finalPipeline).toArray();

        console.log(aggregatedResults);

        // Extract the first (and only) element of the aggregatedResults, which contains your facets and documents
        const aggregationResult = aggregatedResults[0];

        // Extract documents from the 'documents' facet
        const documents = aggregationResult.documents;

        // Extract and transform facets
        const facets = Object.keys(aggregationResult)
            .filter(key => key !== 'documents') // Exclude the 'documents' key to process only facets
            .reduce((acc, key) => {
                // Transform each facet's results for easier consumption
                acc[key] = aggregationResult[key].map(facet => ({
                    _id: facet._id, // Assuming each object has an _id field
                    count: facet.count // Assuming each object has a count field
                }));
                return acc;
            }, {});

        // Assuming totalTasksCount is calculated elsewhere in your code
        res.json({
            results: documents,
            facets: facets,
            pagination: {
                page,
                pageSize,
                pageCount: Math.ceil(totalTasksCount / pageSize),
                totalCount: totalTasksCount,
            }
        });
    } catch (error) {
        console.error('API request failed:', error);
        res.status(500).send('An error occurred while fetching tasks.');
    } finally {
        if (client) {
            await client.close();
        }
    }
});


// Start the server
const PORT = process.env.PORT || 3001;
app.listen(PORT, () => {
    console.log(`Server running on port ${PORT}`);
});