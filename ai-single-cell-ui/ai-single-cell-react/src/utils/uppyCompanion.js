const express = require('express');
const bodyParser = require('body-parser');
const session = require('express-session')
const companion = require('@uppy/companion')

const app = express()

// Companion requires body-parser and express-session middleware.
// You can add it like this if you use those throughout your app.
//
// If you are using something else in your app, you can add these
// middlewares in the same subpath as Companion instead.
app.use(bodyParser.json())
app.use(session({ secret: 'MY_COMPANION_SECRET' }))

const options = {
    providerOptions: {
        drive: {
            key: '848148244825-4sepur31cfcrkc1skh6mt8fiaighv3b2.apps.googleusercontent.com',
            secret: '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',
        },
        onedrive: {
            key: 'my_key',
            secret: 'my_secret',
        },
        dropbox: {
            key: 'my_key',
            secret: 'my_secret',
        }
    },
    server: {
        host: 'localhost:3020',
        protocol: 'http',
        // This MUST match the path you specify in `app.use()` below:
        path: '/companion',
    },
    secret: "MY_COMPANION_SECRET",
    filePath: 'D:/Work/storage/',
}

const { app: companionApp } = companion.app(options)

app.use('/companion', companionApp)