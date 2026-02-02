const express = require('express');
const router = express.Router();
const OpenAI = require('openai');
const { GoogleGenerativeAI } = require("@google/generative-ai");
const multer = require('multer');
const pdf = require('pdf-parse');
const fs = require('fs');

// Configure multer for file uploads (memory storage)
const upload = multer({
    storage: multer.memoryStorage(),
    limits: { fileSize: 10 * 1024 * 1024 } // 10MB limit
});

// Initialize OpenAI client
// Note: It's best practice to use environment variables for API keys
// For now, we'll try to read from process.env, but placeholders are here if needed.
// Users should add OPENAI_API_KEY and GEMINI_API_KEY to their .env file.

router.post('/', upload.single('file'), async (req, res) => {
    let { message, model } = req.body;
    const file = req.file;

    console.log(`[Chat API] Received request for model: ${model}`);

    // Log key presence (checking length to be safe against empty strings)
    console.log(`[Chat API] OpenAI Key present: ${process.env.OPENAI_API_KEY ? 'Yes' : 'No'}`);
    console.log(`[Chat API] Gemini Key present: ${process.env.GEMINI_API_KEY ? 'Yes' : 'No'}`);

    if (!message && !file) {
        return res.status(400).json({ error: 'Message or file is required' });
    }

    // Default message if only file is provided
    if (!message) {
        message = "Please analyze the attached file.";
    }

    try {
        let fileContext = '';
        let imagePart = null;

        if (file) {
            console.log(`[Chat API] Processing file: ${file.originalname} (${file.mimetype})`);

            if (file.mimetype === 'application/pdf') {
                try {
                    const data = await pdf(file.buffer);
                    fileContext = `\n\n[Attached PDF Content: ${file.originalname}]\n${data.text}\n[End of PDF]\n`;
                } catch (e) {
                    console.error("Error parsing PDF:", e);
                    fileContext = `\n\n[Error reading PDF file: ${file.originalname}]`;
                }
            } else if (file.mimetype.startsWith('text/') || file.mimetype === 'application/json' || file.mimetype === 'text/csv' || file.originalname.endsWith('.js') || file.originalname.endsWith('.py')) {
                const text = file.buffer.toString('utf8');
                fileContext = `\n\n[Attached File Content: ${file.originalname}]\n${text}\n[End of File]\n`;
            } else if (file.mimetype.startsWith('image/')) {
                // For images, we process them differently based on the model
                const base64Image = file.buffer.toString('base64');
                imagePart = {
                    mimeType: file.mimetype,
                    data: base64Image
                };
            } else {
                fileContext = `\n\n[Attached File: ${file.originalname} (Unsupported type: ${file.mimetype})]`;
            }
        }

        let reply = '';

        if (model === 'gpt') {
            const openai = new OpenAI({
                apiKey: process.env.OPENAI_API_KEY,
            });

            const messages = [];

            if (imagePart) {
                messages.push({
                    role: "user",
                    content: [
                        { type: "text", text: message },
                        {
                            type: "image_url",
                            image_url: {
                                "url": `data:${imagePart.mimeType};base64,${imagePart.data}`
                            }
                        }
                    ]
                });
            } else {
                messages.push({ role: "user", content: message + fileContext });
            }

            const completion = await openai.chat.completions.create({
                messages: messages,
                model: "gpt-4o", // or gpt-4
            });
            reply = completion.choices[0].message.content;

        } else if (model === 'gemini') {
            const genAI = new GoogleGenerativeAI(process.env.GEMINI_API_KEY);
            // Use gemini-1.5-flash-001 for specific version compatibility
            console.log("Using Gemini Model: gemini-1.5-flash-001");
            const geminiModel = genAI.getGenerativeModel({ model: "gemini-1.5-flash-001" });

            let prompt = message + fileContext;
            let result;

            if (imagePart) {
                // Gemini format for inline data
                const image = {
                    inlineData: {
                        data: imagePart.data,
                        mimeType: imagePart.mimeType,
                    },
                };
                result = await geminiModel.generateContent([prompt, image]);
            } else {
                result = await geminiModel.generateContent(prompt);
            }

            const response = await result.response;
            reply = response.text();

        } else {
            return res.status(400).json({ error: 'Invalid model selection' });
        }

        res.json({ reply });

    } catch (error) {
        console.error('Chat API Error:', error);
        if (error.response) {
            console.error('Chat API Error Response:', error.response.data);
        }

        const status = error.status || 500;
        // Return more specific error message if possible
        const errorMessage = error.message || 'Failed to fetch response from AI provider';
        res.status(status).json({ error: 'AI Provider Error', details: errorMessage });
    }
});

module.exports = router;
