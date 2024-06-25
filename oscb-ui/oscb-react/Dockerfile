FROM node:14
# FROM stormgod/oscb_react

# Create app directory
WORKDIR /opt/src/app
# Install app dependencies
COPY package*.json ./

RUN npm install
# RUN npm install @octokit/rest
# RUN npm install @rjsf/utils
# RUN npm install @rjsf/validator-ajv8
# Copy app source code
COPY . .

#Expose port and start application
EXPOSE 3000
CMD ["npm", "start"]
