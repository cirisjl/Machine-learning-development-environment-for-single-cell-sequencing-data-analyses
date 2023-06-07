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

console.log('HOSTURL: ' + process.env.HOST_URL);
const app = express();
app.use(cors({
    origin: [`http://${process.env.HOST_URL}:3000`, `http://${hostIp}:3000`],
    credentials: true
}));
app.use(bodyParser.json());
app.use(cookieParser());

const dbConfig = JSON.parse(fs.readFileSync('./configs/dbconfigs.json'));
const storageConfig = JSON.parse(fs.readFileSync('./configs/storageConfig.json'));
const { storageDir, storageAllowance, intermediateStorage } = storageConfig;

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
                    const jwtToken = jwt.sign({ username, password }, 'secret', { expiresIn: '1h' });

                    // the cookie will be set with the name "jwtToken" and the value of the token
                    // the "httpOnly" and "secure" options help prevent XSS and cookie theft
                    // the "secure" option is only set if the app is running in production mode
                    // set the cookie with the JWT token on the response object
                    res.cookie("jwtToken", jwtToken, {
                        //httpOnly: true,
                        maxAge: 60 * 60 * 1000,
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
            const jwtToken = jwt.sign({ username, password }, 'secret', { expiresIn: '1h' });

            // the cookie will be set with the name "jwtToken" and the value of the token
            // the "httpOnly" and "secure" options help prevent XSS and cookie theft
            // the "secure" option is only set if the app is running in production mode
            // set the cookie with the JWT token on the response object
            res.cookie("jwtToken", jwtToken, {
                //httpOnly: true,
                maxAge: 60 * 60 * 1000,
                path: "/"
                //secure: process.env.NODE_ENV === "production",
            });

            console.log(req.cookies);
            res.json({ status: 200, message: 'Logged in successfully', jwtToken });
        });
    });
});

// Route to handle protected resource
app.get('/protected', verifyToken, (req, res) => {
    jwt.verify(req.token, 'secret', (err, authData) => {
        if (err) {
            res.sendStatus(403);
        } else {
            res.json({ message: 'You have access to the protected resource', authData });
        }
    });
});

app.post('/createDataset', (req, res) => {

    const { title, n_cells, reference, summary, authToken, files } = req.body;
    const username = getUserFromToken(authToken);

    pool.getConnection(function (err, connection) {
        if (err) throw err;

        connection.beginTransaction(function (err) {
            if (err) throw err;

            // Run SELECT command
            connection.query('SELECT user_id FROM users WHERE username = ? LIMIT 1', [username], function (err, userRows) {
                if (err) {
                    connection.rollback(function () {
                        throw err;
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
                                throw err;
                            });
                        }
                    } else {
                        const datasetId = datasetResult.insertId;

                        for (const file of files) {
                            connection.query('INSERT INTO file (file_loc, dataset_id) VALUES (?, ?)', [file, datasetId]);
                        }

                        // Commit transaction
                        connection.commit(function (err) {
                            if (err) {
                                connection.rollback(function () {
                                    throw err;
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
});

app.put('/updateDataset', async (req, res) => {

    const { title, n_cells, reference, summary, authToken, files, currentFileList } = req.body;
    const username = getUserFromToken(authToken);

    const insertList = files.filter(item => !currentFileList.includes(item));
    const deleteList = currentFileList.filter(item => !files.includes(item));

    pool.getConnection(function (err, connection) {
        if (err) throw err;

        connection.beginTransaction(function (err) {
            if (err) { console.log('Idi error: ' + err); throw err; }

            // Run SELECT command
            connection.query('SELECT user_id FROM users WHERE username = ? LIMIT 1', [username], function (err, userRows) {
                if (err) {
                    connection.rollback(function () {
                        throw err;
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

                connection.query('SELECT dataset_id FROM dataset WHERE user_id = ? and title = ? LIMIT 1', [userId, title], function (err, datasetRows) {
                    if (err) {
                        connection.rollback(function () {
                            throw err;
                        });
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
                                throw err;
                            });
                        }
                        for (const file of insertList) {
                            connection.query('INSERT INTO file (file_loc, dataset_id) VALUES (?, ?)', [file, datasetId]);
                        }

                        for (const file of deleteList) {
                            connection.query('DELETE FROM file WHERE file_loc=? AND dataset_id=?', [file, datasetId]);
                        }

                        // Commit transaction
                        connection.commit(function (err) {
                            if (err) {
                                connection.rollback(function () {
                                    throw err;
                                });
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
        if (err) throw err;

        connection.beginTransaction(function (err) {
            if (err) { throw err; }

            // Run SELECT command
            connection.query('SELECT user_id FROM users WHERE username = ? LIMIT 1', [username], function (err, userRows) {
                if (err) {
                    connection.rollback(function () {
                        throw err;
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

                connection.query('SELECT dataset_id FROM dataset WHERE user_id = ? and title = ? LIMIT 1', [userId, dataset], function (err, datasetRows) {
                    if (err) {
                        connection.rollback(function () {
                            throw err;
                        });
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

                            connection.rollback(function () {
                                throw err;
                            });
                        }
                        connection.query('DELETE FROM dataset where dataset_id=?', [datasetId]);

                        // Commit transaction
                        connection.commit(function (err) {
                            if (err) {
                                connection.rollback(function () {
                                    throw err;
                                });
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

    pool.query(`UPDATE aisinglecell.file f JOIN aisinglecell.dataset d ON f.dataset_id = d.dataset_id JOIN aisinglecell.users u ON d.user_id = u.user_id SET f.file_loc = CONCAT('${newName}', SUBSTRING(f.file_loc, LENGTH('${oldName}') + 1)) WHERE u.username = '${uname}' AND f.file_loc LIKE '${oldName}%';`, [newName, oldName, uname, oldName + '%'], (err, results) => {

        console.log(results);
        if (err) {
            console.error(err);
            res.json({ status: 500, message: 'Internal Server Error' });
            return;
        }
    })

    fs.rename(`${storageDir}${uname}/${oldName}`, `${storageDir}${uname}/${newName}`, (err) => {
        if (err) {
            console.error(err);
        } else {
            console.log('File renamed successfully!');
        }
    });
    return res.status(200).jsonp('Ok');
});

app.post('/download', async (req, res) => {
    const { fileList } = req.body;
    const { authToken } = req.query;
    console.log('Entered download function');

    const username = getUserFromToken(authToken);
    console.log('Username is: ' + username);
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
            const filePath = path.join(storageDir, username, item);
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
    const username = getUserFromToken(authToken);
    let filePath = '';
    console.log('Entered download function');

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
    console.log('Entered download function');

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
                const [rows, fields] = await pool.promise().execute(`SELECT f.file_loc FROM aisinglecell.file f JOIN aisinglecell.dataset d ON f.dataset_id = d.dataset_id JOIN aisinglecell.users u ON d.user_id = u.user_id WHERE u.username = '${uname}' AND f.file_loc like '${file.replace("'", "''")}%';`);
                console.log('Count: ' + rows.length);
                if (rows.length > 0) {
                    res.status(401).json({ message: 'File(s) being used by datasets.' });
                    return;
                } else {
                    const filePath = `${storageDir}${uname}/${file}`;
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

        if (usingFor === 'resultFiles')
            var directoryPath = path.join(intermediateStorage + "/" + dirPath + "/");
        else {
            var directoryPath = path.join(storageDir + uid + "/" + dirPath + "/");
            if (subdir != undefined)
                directoryPath = path.join(storageDir + uid + "/", subdir);
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
    let { uploadDir, authToken } = req.query;
    let username = getUserFromToken(authToken);

    let destDir = `./storage/${username}/${uploadDir}`; // Replace with your storage directory
    let tempDir = './uploads'; // Replace with a temporary directory for uploads

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
    const folderPath = `${storageDir}/${username}/${pwd}/${folderName}`;
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
        if (err) throw err;

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
                if (err) throw err;

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
                            direc
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
    const { authToken, workflow, taskId } = req.body;
    const username = getUserFromToken(authToken);

    pool.getConnection(function (err, connection) {
        if (err) throw err;

        connection.beginTransaction(function (err) {
            if (err) throw err;

            connection.query('SELECT user_id FROM users WHERE username = ? LIMIT 1', [username], function (err, userRows) {
                if (err) {
                    connection.rollback(function () {
                        throw err;
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

                const date = new Date();
                const timestamp = Date.UTC(date.getUTCFullYear(), date.getUTCMonth(), date.getUTCDate(), date.getUTCHours(), date.getUTCMinutes(), date.getUTCSeconds(), date.getUTCMilliseconds());
                connection.query('INSERT INTO task (user_id, dataset_id, task_id, workflow, created_datetime) VALUES (?, ?, ?, ?, ?)', [userId, 31, taskId, workflow, timestamp], function (err, taskResult) {
                    if (err) {
                        connection.rollback(function () {
                            throw err;
                        });
                    } else {
                        const taskId = taskResult.insertId;

                        // Commit transaction
                        connection.commit(function (err) {
                            if (err) {
                                connection.rollback(function () {
                                    throw err;
                                });
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
        if (err) throw err;

        connection.beginTransaction(function (err) {
            if (err) throw err;

            const date = new Date();
            const timestamp = Date.UTC(date.getUTCFullYear(), date.getUTCMonth(), date.getUTCDate(), date.getUTCHours(), date.getUTCMinutes(), date.getUTCSeconds(), date.getUTCMilliseconds());
            connection.query('UPDATE task SET status = ?, finish_datetime = ? WHERE task_id IN (?)', [status, timestamp, taskIdsArr], function (err, taskResult) {
                if (err) {
                    connection.rollback(function () {
                        throw err;
                    });
                } else {
                    // Commit transaction
                    connection.commit(function (err) {
                        if (err) {
                            connection.rollback(function () {
                                throw err;
                            });
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
        if (err) throw err;

        connection.beginTransaction(function (err) {
            if (err) throw err;

            connection.query('SELECT user_id FROM users WHERE username = ? LIMIT 1', [username], function (err, userRows) {
                if (err) {
                    connection.rollback(function () {
                        throw err;
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

                connection.query('SELECT task_id, workflow, dataset_id, status, created_datetime, finish_datetime FROM task WHERE user_id = ?', [userId], function (err, rows) {
                    if (err) {
                        connection.rollback(function () {
                            throw err;
                        });
                    } else {
                        connection.release();
                        res.json(rows);
                    }
                });
            });
        });
    });
});
// Start the server
const PORT = process.env.PORT || 3001;
app.listen(PORT, () => {
    console.log(`Server running on port ${PORT}`);
});