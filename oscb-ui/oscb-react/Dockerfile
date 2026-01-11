# FROM node:22.21.1
FROM stormgod/oscb-react:0.1

# Create app directory
WORKDIR /opt/src/app
# Install app dependencies
# COPY package*.json ./

# RUN npm config set legacy-peer-deps true
# RUN npm install
# RUN npm install @octokit/rest
# RUN npm install @rjsf/core @rjsf/utils @rjsf/validator-ajv8 --save

# RUN npm install @mui/icons-material
# RUN npm install @mui/material @emotion/react @emotion/styled


# Copy app source code
COPY . .

#Expose port and start application
EXPOSE 3000
CMD ["npm", "start"]
