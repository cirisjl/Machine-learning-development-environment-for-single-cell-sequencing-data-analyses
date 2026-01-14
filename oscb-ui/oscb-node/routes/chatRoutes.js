// Polyfill for fetch which is required for openai in node < 18
const fetch = require('node-fetch');
const FormData = require('form-data');
if (!globalThis.fetch) {
    globalThis.fetch = fetch;
    globalThis.Headers = fetch.Headers;
    globalThis.Request = fetch.Request;
    globalThis.Response = fetch.Response;
    globalThis.FormData = FormData;
    if (!globalThis.Blob) {
        globalThis.Blob = class Blob {
            constructor(content, options) {
                this.content = content;
                this.options = options;
            }
        };
    }
    globalThis.File = class File extends globalThis.Blob {
        constructor(parts, filename, options) {
            super(parts, options);
            this.name = filename;
        }
    }
}

const express = require('express');
const router = express.Router();
const OpenAI = require('openai');
const { GoogleGenerativeAI } = require("@google/generative-ai");

// Initialize OpenAI client
// Note: It's best practice to use environment variables for API keys
// For now, we'll try to read from process.env, but placeholders are here if needed.
// Users should add OPENAI_API_KEY and GEMINI_API_KEY to their .env file.

router.post('/', async (req, res) => {
    const { message, model } = req.body;
    console.log(`[Chat API] Received request for model: ${model}`);
    console.log(`[Chat API] OpenAI Key exists: ${!!process.env.OPENAI_API_KEY}`);
    console.log(`[Chat API] Gemini Key exists: ${!!process.env.GEMINI_API_KEY}`);

    if (!message) {
        return res.status(400).json({ error: 'Message is required' });
    }

    try {
        let reply = '';

        if (model === 'gpt') {
            const openai = new OpenAI({
                apiKey: process.env.OPENAI_API_KEY,
            });

            const completion = await openai.chat.completions.create({
                messages: [{ role: "user", content: message }],
                model: "gpt-3.5-turbo", // or gpt-4
            });
            reply = completion.choices[0].message.content;

        } else if (model === 'gemini') {
            const genAI = new GoogleGenerativeAI(process.env.GEMINI_API_KEY);
            const geminiModel = genAI.getGenerativeModel({ model: "gemini-2.0-flash" });

            const result = await geminiModel.generateContent(message);
            const response = await result.response;
            reply = response.text();

        } else {
            return res.status(400).json({ error: 'Invalid model selection' });
        }

        res.json({ reply });

    } catch (error) {
        console.error('Chat API Error:', error);
        res.status(500).json({ error: 'Failed to fetch response from AI provider', details: error.message });
    }
});

module.exports = router;
