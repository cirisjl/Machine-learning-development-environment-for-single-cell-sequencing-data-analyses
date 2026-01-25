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
    
    // Log key presence (checking length to be safe against empty strings)
    console.log(`[Chat API] OpenAI Key present: ${process.env.OPENAI_API_KEY ? 'Yes' : 'No'}`);
    console.log(`[Chat API] Gemini Key present: ${process.env.GEMINI_API_KEY ? 'Yes' : 'No'}`);

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
                model: "gpt-4o", // or gpt-4
            });
            reply = completion.choices[0].message.content;

        } else if (model === 'gemini') {
            const genAI = new GoogleGenerativeAI(process.env.GEMINI_API_KEY);
            const geminiModel = genAI.getGenerativeModel({ model: "gemini-flash-latest" });

            const result = await geminiModel.generateContent(message);
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
