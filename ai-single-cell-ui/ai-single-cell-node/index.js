const express = require('express');
const fs = require('fs-extra');
const mysql = require('mysql2/promise');
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

let mongoDBConfig; // Declare mongoDBConfig outside the try-catch block
let mongoUrl, dbName, optionsCollectionName, datasetCollectionName;

try {
    mongoDBConfig = JSON.parse(fs.readFileSync('./configs/mongoDB.json')); // Import the MongoDB connection configuration
} catch (error) {
    console.error('Error reading MongoDB configuration:', error);
}

// Check if mongoDBConfig is defined before using it
if (mongoDBConfig) {
    ({ mongoUrl, dbName, optionsCollectionName, datasetCollectionName } = mongoDBConfig);
    const { MongoClient, ObjectId } = require('mongodb');
    // Now you can use MongoClient and ObjectId here
} else {
    console.error('MongoDB configuration not available. Ensure ./configs/mongoDB.json is properly configured.');
}

console.log('HOSTURL: ' + process.env.HOST_URL);
const app = express();
app.use(cors({
    origin: [`http://${process.env.HOST_URL}:3000`, `http://${hostIp}:3000`, 'http://node-0.jiangl0-160204.biomizzou-pg0.clemson.cloudlab.us:3000','http://node-0.jiangl0-161295.biomizzou-pg0.clemson.cloudlab.us:3000','http://node-0.ai-single-cell.biomizzou-pg0.clemson.cloudlab.us:3000'],
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
        const decoded = jwt.verify(token, 'secret');

        console.log('token: ' + decoded.username);
        if (!decoded.username) {
            return 'Unauthorized';
        }

        return decoded.username;
    } catch (err) {
        console.log('eee ' + err)
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
    const targetFolderName = index === 1 ?  path.join(destinationDir, folderName) : path.join(destinationDir, `${folderName}(${index})`);
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

        // if(dirName) {
        //     destinationDir = path.join(destinationDir, dirName);
        // } 
        console.log("logger to debug the source and destination directories");
        console.log("source" + sourceDir);
        console.log(destinationDir);

        // Ensure the destination directory exists before copying files
    //   await createDirectoryIfNotExists(destinationDir);
    //   const files = await fs.readdir(sourceDir);
  
      for (let file of files) {
        const sourceFilePath = path.join(sourceDir, file);
        let destinationFilePath = "";
        if(fromPublic) {
            file = file.replace(/^\/?publicDatasets\//, '/'); // Remove "PUBLIC_DATASETS" from the start
            destinationFilePath = path.join(destinationDir, file);
        } else {
            destinationFilePath = path.join(destinationDir, file);
        }

        const sourceFileDir = path.dirname(sourceFilePath);
        const destinationFileDir = path.dirname(destinationFilePath);
  
        // Ensure the destination directory exists before copying files
        await createDirectoryIfNotExists(destinationFileDir);
 
        console.log("final paths");
        console.log("source paths" + sourceFilePath);
        console.log("destination paths" + destinationFilePath);

        // Perform the actual file copy
        await fs.copyFile(sourceFilePath, destinationFilePath);
      }
    } catch (error) {
      console.error('Error copying files:', error);
    }
  };

  app.post('/api/signup', async (req, res) => {
    const { username, email, password } = req.body;

    try {
        // Hash the password using bcrypt
        const hash = await bcrypt.hash(password, 10);

        // Insert the user into the database
        await pool.query('INSERT INTO users (username, email, password_hash) VALUES (?, ?, ?)', [username, email, hash]);

        // Create directory for user storage
        if (!fs.existsSync(storageDir + username)) {
            await fs.promises.mkdir(storageDir + username);
        }

        // Create JWT token and send it back to the client
        const jwtToken = jwt.sign({ username }, 'secret' , { expiresIn: '1h' });

        res.cookie("jwtToken", jwtToken, {
            httpOnly: true,
            maxAge: 60 * 60 * 1000,
            path: "/",
            // secure: process.env.NODE_ENV === "production",
        });

        res.json({ status: 200, message: 'User account created successfully' });
    } catch (error) {
        console.error('Signup error:', error);

        if (error.code === 'ER_DUP_ENTRY') {
            res.status(409).json({ status: 409, message: 'Username or email already exists' });
        } else {
            res.status(500).json({ status: 500, message: 'Internal Server Error' });
        }
    }
});


app.post('/api/login', async (req, res) => {
    const { username, password } = req.body;

    try {
        // Query the database for the user
        const [results] = await pool.query('SELECT * FROM users WHERE username = ?', [username]);

        if (results.length === 0) {
            return res.status(401).json({ message: 'Invalid credentials' });
        }

        const user = results[0];

        // Compare the provided password with the stored hash
        const isMatch = await bcrypt.compare(password, user.password_hash);
        if (!isMatch) {
            return res.status(401).json({ message: 'Invalid credentials' });
        }

        // Create JWT token
        const jwtToken = jwt.sign({ username }, 'secret', { expiresIn: '1h' });

        // Set the cookie with the JWT token
        res.cookie("jwtToken", jwtToken, {
            // httpOnly: true,
            maxAge: 60 * 60 * 1000, // 1 hour
            path: "/",
            // secure: process.env.NODE_ENV === "production",
        });

        res.json({ status: 200, message: 'Logged in successfully', jwtToken });
    } catch (error) {
        console.error('Login error:', error);
        res.status(500).json({ message: 'Internal Server Error' });
    }
});

app.get('/protected', verifyToken, async (req, res) => {
    try {
        const authData = jwt.verify(req.token, 'secret');

        if (!authData.username) {
            return res.status(403).json({ message: 'Access Denied' });
        }

        const [results] = await pool.query('SELECT isAdmin FROM users WHERE username = ?', authData.username);

        if (results.length === 0) {
            return res.status(401).json({ message: 'Invalid credentials' });
        }

        authData.isAdmin = results[0].isAdmin === 1;
        res.json({ message: 'You have access to the protected resource', authData });

    } catch (error) {
        if (error.name === 'JsonWebTokenError') {
            res.status(403).json({ message: 'Invalid Token' });
        } else {
            console.error('Protected route error:', error);
            res.status(500).json({ message: 'Internal Server Error' });
        }
    }
});

app.post('/createDataset', async (req, res) => {

    const { title, n_cells, reference, summary, authToken, files, makeItpublic } = req.body;
    const username = getUserFromToken(authToken);

    let filesFromPublic = false;

    console.log("Logger to debug makeit public flag : " + makeItpublic);

    // Logic to Copy files from public storage to user private storage if it is a public Dataset.
    for (const file of files) {
        console.log("inside for loop")
        console.log(file);
        if(file.startsWith("publicDataset") || file.startsWith("/publicDatasets")) {
            console.log("inside if loop for my check");
            filesFromPublic = true;
            break;
        }
    }
    console.log("value of filesFromPublic::: " + filesFromPublic);

    if(filesFromPublic) {
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
                            if(filesFromPublic) {
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

    if(makeItpublic) {
        try {
            let dirName = "";
            const fromPublic = false;
            if (files.length > 0) {
                dirName = path.dirname(files[0])
            } 

            let userPrivateStorageDir = storageDir + username // Change this to the user's private storage path

            console.log("logger to see the userPrivateStorageDir" + userPrivateStorageDir);

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

    let filesFromPublic = files.some(file => file.startsWith("publicDatasets") || file.startsWith("/publicDatasets"));

    if (filesFromPublic) {
        try {
            let dirName = files.length > 0 ? path.dirname(files[0]) : "";
            let userPrivateStorageDir = storageDir + username;
            await copyFiles("/usr/src/app/storage/", userPrivateStorageDir, dirName, files, filesFromPublic);
        } catch (err) {
            console.error('Error copying files:', err);
            return res.status(500).send('Error in file operation');
        }
    }

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

            connection.query('SELECT user_id FROM users WHERE username = ? LIMIT 1', [username], function (err, userRows) {
                if (err) {
                    console.error('Error in SELECT query:', err);
                    connection.rollback(function () {
                        connection.release();
                    });
                    return res.status(500).send('Database query error');
                }

                const userId = userRows[0]?.user_id;
                if (!userId) {
                    res.status(400).send('User not found');
                    connection.rollback(function () {
                        connection.release();
                    });
                    return;
                }

                connection.query('SELECT dataset_id FROM dataset WHERE user_id = ? and title = ? LIMIT 1', [userId, title], function (err, datasetRows) {
                    if (err) {
                        console.error('Error in SELECT query for dataset:', err);
                        connection.rollback(function () {
                            connection.release();
                        });
                        return res.status(500).send('Dataset query error');
                    }

                    const datasetId = datasetRows[0]?.dataset_id;
                    if (!datasetId) {
                        res.status(400).send('Dataset not found');
                        connection.rollback(function () {
                            connection.release();
                        });
                        return;
                    }

                    connection.query('UPDATE dataset SET n_cells=?, reference=?, summary=? WHERE dataset_id=?', [n_cells, reference, summary, datasetId], async function (err, datasetResult) {
                        if (err) {
                            console.error('Error in UPDATE query:', err);
                            connection.rollback(function () {
                                connection.release();
                            });
                            return res.status(500).send('Dataset update error');
                        }

                        // Insert new files
                        try {
                            for (let file of insertList) {
                                file = filesFromPublic ? file.replace(/^\/?publicDatasets\//, '/') : file;
                                await connection.query('INSERT INTO file (file_loc, dataset_id) VALUES (?, ?)', [file, datasetId]);
                            }

                            // Delete removed files
                            for (const file of deleteList) {
                                await connection.query('DELETE FROM file WHERE file_loc=? AND dataset_id=?', [file, datasetId]);
                            }
                        } catch (err) {
                            console.error('Error in file insert/delete:', err);
                            connection.rollback(function () {
                                connection.release();
                            });
                            return res.status(500).send('File insert/delete error');
                        }

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

            connection.query('SELECT user_id FROM users WHERE username = ? LIMIT 1', [username], function (err, userRows) {
                if (err) {
                    console.error('Error in SELECT query:', err);
                    connection.rollback(function () {
                        connection.release();
                    });
                    return res.status(500).send('Database query error');
                }

                const userId = userRows[0]?.user_id;
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

                    const datasetId = datasetRows[0]?.dataset_id;
                    if (!datasetId) {
                        res.status(400).send('Dataset not found');
                        connection.rollback(function () {
                            connection.release();
                        });
                        return;
                    }

                    connection.query('DELETE FROM file WHERE dataset_id=?', [datasetId], function (err) {
                        if (err) {
                            console.error('Error deleting files:', err);
                            connection.rollback(function () {
                                connection.release();
                            });
                            return res.status(500).send('Error deleting files');
                        }

                        connection.query('DELETE FROM dataset WHERE dataset_id=?', [datasetId], function (err) {
                            if (err) {
                                console.error('Error deleting dataset:', err);
                                connection.rollback(function () {
                                    connection.release();
                                });
                                return res.status(500).send('Error deleting dataset');
                            }

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
});

app.post('/renameFile', async (req, res) => {
    let { oldName, newName, authToken } = req.query;

    oldName = oldName.replace('//', '/');
    newName = newName.replace('//', '/');

    const uname = getUserFromToken(authToken);
    if (uname == 'Unauthorized') {
        return res.status(403).jsonp('Unauthorized');
    }

    try {
        // Avoiding SQL injection by using prepared statements
        const query = `SELECT f.file_loc FROM aisinglecell.file f JOIN aisinglecell.dataset d ON f.dataset_id = d.dataset_id JOIN aisinglecell.users u ON d.user_id = u.user_id WHERE u.username = ? AND f.file_loc LIKE ?;`;
        const [results] = await pool.query(query, [uname, `${oldName}%`]);

        if (results.length > 0) {
            return res.status(409).json({ message: 'Directory already exists' });
        } else {
            const oldPath = oldName.includes("publicDatasets") ? `/usr/src/app/storage/${oldName}` : `${storageDir}${uname}/${oldName}`;
            const newPath = newName.includes("publicDatasets") ? `/usr/src/app/storage/${newName}` : `${storageDir}${uname}/${newName}`;

            await fs.promises.rename(oldPath, newPath);
            console.log('File renamed successfully!');
            res.status(200).jsonp('Ok');
        }
    } catch (err) {
        console.error(err);
        res.status(500).json({ message: 'Internal Server Error' });
    }
});


app.post('/download', async (req, res) => {
    const { fileList } = req.body;
    const { authToken, pwd } = req.query;

    const username = getUserFromToken(authToken);
    if (username === 'Unauthorized') {
        return res.status(403).json({ message: 'Unauthorized' });
    }

    if (!fileList || !Array.isArray(fileList) || fileList.length === 0) {
        return res.status(400).json({ message: 'Invalid request' });
    }

    try {
        const zipName = 'files.zip';
        const output = fs.createWriteStream(zipName);
        const archive = archiver('zip'); // Sets the compression level

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
            if(pwd && pwd.includes("publicDatasets")) {
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

        archive.on('error', (err) => {
            console.error('Archive error:', err);
            res.status(500).json({ message: 'Error creating archive' });
        });

    } catch (err) {
        console.error('Error in download route:', err);
        res.status(500).json({ message: 'Internal Server Error' });
    }
});

app.get('/download', async (req, res) => {
    const { fileUrl, authToken, forResultFile } = req.query; // forResultFile seems unused
    const { pwd } = req.query;
    const username = getUserFromToken(authToken);

    if (!fileUrl || username === 'Unauthorized') {
        return res.status(400).jsonp('Invalid request or Unauthorized');
    }

    let filePath = pwd && pwd.includes("publicDatasets") ? path.join(storageDir, fileUrl) : path.join(storageDir, username, fileUrl);
    console.log('file: ' + filePath);

    try {
        const fileStat = await fs.promises.stat(filePath);

        if (fileStat.isFile()) {
            // Download file
            const filename = path.basename(fileUrl);
            const mimetype = mime.getType(filePath, { legacy: true }) || 'application/octet-stream';

            res.setHeader('Content-Disposition', `attachment; filename="${filename}"`);
            res.setHeader('Content-Type', mimetype);

            const filestream = fs.createReadStream(filePath);
            filestream.pipe(res);

            filestream.on('error', (err) => {
                console.error('File stream error:', err);
                res.status(500).jsonp('Internal Server Error');
            });
        } else if (fileStat.isDirectory()) {
            // Download folder as zip
            const folderName = path.basename(filePath);
            const archive = archiver('zip', { zlib: { level: 9 } });

            archive.on('error', (err) => {
                console.error('Archive error:', err);
                res.status(500).jsonp('Internal Server Error');
            });

            res.setHeader('Content-Disposition', `attachment; filename="${folderName}.zip"`);
            res.setHeader('Content-Type', 'application/zip');

            archive.directory(filePath, folderName);
            archive.pipe(res);
            archive.finalize();
        } else {
            return res.status(400).jsonp('Invalid request');
        }
    } catch (err) {
        console.error('Error in download route:', err);
        res.status(500).jsonp('Internal Server Error');
    }
});

app.get('/fetchPreview', async (req, res) => {
    const { fileUrl, authToken, forResultFile } = req.query;
    const username = getUserFromToken(authToken);

    if (!fileUrl || username === 'Unauthorized') {
        return res.status(400).jsonp('Invalid request or Unauthorized');
    }

    let filePath = forResultFile ? `${intermediateStorage}/${fileUrl}` : `${storageDir}/${username}/${fileUrl}`;
    console.log('file: ' + filePath);

    try {
        const fileStat = await fs.promises.stat(filePath);

        if (fileStat.isFile()) {
            // Read first 20 lines of the file
            const fileStream = fs.createReadStream(filePath, { encoding: 'utf8' });
            let lines = '';
            let lineCount = 0;

            fileStream.on('data', (chunk) => {
                lines += chunk;
                lineCount = lines.split('\n').length;

                if (lineCount >= 20) {
                    fileStream.destroy();
                    lines = lines.split('\n').slice(0, 20).join('\n');
                    res.status(200).send(lines);
                }
            });

            fileStream.on('end', () => {
                if (lineCount < 20) {
                    res.status(200).send(lines);
                }
            });

            fileStream.on('error', (error) => {
                console.error('Error reading file:', error);
                res.status(500).jsonp('Error reading file');
            });
        } else {
            res.status(400).jsonp('File is not accessible');
        }
    } catch (error) {
        console.error('Error in fetchPreview route:', error);
        res.status(500).jsonp('Internal Server Error');
    }
});


app.delete('/deleteFiles', async (req, res) => {
    const { fileList } = req.body;
    const { authToken } = req.query;
    const { pwd } = req.query
    console.log('Entered delete function');
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
                const [rows, fields] = await pool.execute(`SELECT f.file_loc FROM aisinglecell.file f JOIN aisinglecell.dataset d ON f.dataset_id = d.dataset_id JOIN aisinglecell.users u ON d.user_id = u.user_id WHERE u.username = '${uname}' AND f.file_loc like '${file.replace("'", "''")}%';`);
                console.log('Count: ' + rows.length);
                if (rows.length > 0) {
                    res.status(401).json({ message: 'File(s) being used by datasets.' });
                    return;
                } else {
                    let filePath = ""
                    if(pwd.includes("publicDatasets")) {
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
        console.log("Debug point to check the value of subdir");
        console.log(dirPath);
        var directoryPath = ""

        
        var directoryPath = path.join(storageDir + uid + "/" + dirPath + "/");
        
        if (subdir != undefined)
            directoryPath = path.join(storageDir + uid + "/", subdir);

        if(dirPath == "publicDatasets") {
            directoryPath = publicStorage;
        }

        if(dirPath.includes("publicDatasets/")) {
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
    let { uploadDir, authToken ,publicDatasetFlag} = req.query;
    let username = getUserFromToken(authToken);

    let destDir = publicDatasetFlag === "true" ? "./storage/" + uploadDir : "./storage/" + username + uploadDir ;
    console.log("publicdatasetFlag value debug point:::: " + publicDatasetFlag);
    console.log("Inside else destDir:: " + destDir);

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

app.post('/createNewFolder', async (req, res) => {
    const { pwd, folderName, authToken } = req.query;
    const username = getUserFromToken(authToken);

    if (username === 'Unauthorized') {
        return res.status(403).json('Unauthorized');
    }

    let folderPath = '';
    if (pwd.includes("publicDatasets")) {
        folderPath = path.join("/usr/src/app/storage", pwd, folderName);
    } else {
        folderPath = path.join(storageDir, username, pwd, folderName);
    }

    try {
        if (fs.existsSync(folderPath)) {
            return res.status(400).jsonp('Folder already exists');
        }

        await fs.promises.mkdir(folderPath, { recursive: true });
        res.status(201).jsonp('Folder created');
    } catch (err) {
        console.error('Error creating folder:', err);
        res.status(500).jsonp('Error creating folder: ' + err.message);
    }
});

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

app.get('/preview/datasets', (req, res) => {
    const { authToken } = req.query;
    const username = getUserFromToken(authToken);

    if (username === "Unauthorized") {
        return res.status(403).json('Unauthorized');
    }

    const userQuery = 'SELECT user_id FROM users WHERE username = ?';

    pool.query(userQuery, [username], (err, userResult) => {
        if (err) {
            console.error('Database query error:', err);
            return res.status(500).json({ message: 'Internal Server Error' });
        }

        if (userResult.length === 0) {
            return res.status(404).json({ message: `User '${username}' not found` });
        }

        const userID = userResult[0].user_id;
        const datasetsQuery = `
            SELECT dataset.dataset_id, dataset.title, dataset.n_cells, dataset.reference, dataset.summary, file.file_id, file.file_loc, SUBSTRING_INDEX(SUBSTRING_INDEX(file.file_loc, '/', 2), '/', -1) AS direc
            FROM dataset
            JOIN file ON dataset.dataset_id = file.dataset_id
            WHERE dataset.user_id = ?
        `;

        pool.query(datasetsQuery, [userID], (err, datasetsResult) => {
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

    if (username === "Unauthorized") {
        return res.status(403).json({ message: 'Unauthorized' });
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

            connection.query('SELECT user_id FROM users WHERE username = ? LIMIT 1', [username], function (err, userRows) {
                if (err) {
                    console.error('Error in SELECT query:', err);
                    connection.rollback(function () {
                        connection.release();
                    });
                    return res.status(500).json({ message: 'Database query error' });
                }

                const userId = userRows[0]?.user_id;
                if (!userId) {
                    res.status(400).json({ message: 'User not found' });
                    connection.rollback(function () {
                        connection.release();
                    });
                    return;
                }

                const timestamp = new Date().toISOString();
                connection.query('INSERT INTO task (task_title, task_id, user_id, tool, results_path, created_datetime) VALUES (?, ?, ?, ?, ?, ?)', [taskTitle, taskId, userId, method, outputPath, timestamp], function (err, taskResult) {
                    if (err) {
                        console.error('Error in INSERT query:', err);
                        connection.rollback(function () {
                            connection.release();
                        });
                        return res.status(500).json({ message: 'Database insertion error' });
                    }

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
                        res.status(201).json({ message: 'Task Created.' });
                    });
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

    if (username === "Unauthorized") {
        return res.status(403).json({ message: 'Unauthorized' });
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

            connection.query('SELECT user_id FROM users WHERE username = ? LIMIT 1', [username], function (err, userRows) {
                if (err) {
                    console.error('Database query error:', err);
                    return res.status(500).json({ message: 'Internal Server Error' });
                }
        
                if (userRows.length === 0) {
                    return res.status(404).json({ message: 'User not found' });
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
    try {
      const client = new MongoClient(mongoUrl, { useUnifiedTopology: true });
  
      // Connect to the MongoDB server
      await client.connect();
  
      const db = client.db(dbName);
      const collection = db.collection(datasetCollectionName);
  
      const formData = req.body; // This assumes you have the necessary middleware to parse JSON in the request body
  
      // Check if a document with the provided Id already exists
      const existingDocument = await collection.findOne({ Id: formData.Id });
  
      if (existingDocument) {
        // Document with the provided Id already exists
        console.log('Document with Id already exists:', formData.Id);
        res.status(400).json({ error: 'Document with the provided Id already exists' });
      } else {
        // Document with the provided Id does not exist, proceed with insertion
        collection.insertOne(formData, (err) => {
          if (err) {
            console.error('Error inserting form data into MongoDB:', err);
            res.status(500).json({ error: 'Error submitting form data' });
          } else {
            console.log('Form data submitted successfully');
            res.status(200).json({ message: 'Form data submitted successfully' });
          }
  
          // Close the MongoDB connection here after the operation is complete
          client.close();
        });
      }
    } catch (err) {
      console.error('Error:', err);
      res.status(500).json({ error: 'Internal Server Error' });
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
        const matchStage = isAdmin=== 'true' ? {} : { username: username };

        console.log(isAdmin);
        console.log(matchStage);

        // Aggregation pipeline stages
        const pipeline = [
            { $match: matchStage },
            {
                $group: {
                    _id: '$field',
                    options: { $addToSet: { _id: '$_id', name: '$name', username: '$username', abbreviation:'$abbreviation' } },
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
    if(username) {
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


// Start the server
const PORT = process.env.PORT || 3001;
app.listen(PORT, () => {
    console.log(`Server running on port ${PORT}`);
});