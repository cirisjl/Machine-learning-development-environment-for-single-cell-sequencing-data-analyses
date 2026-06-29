# FROM node:22.21.1
FROM stormgod/oscb-react:0.1.1

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
RUN npm install react-app-rewired framer-motion react-markdown remark-math rehype-katex katex

# The dependency lockfile belongs to the prebuilt base image, so this derived
# image cannot persist a Browserslist database refresh. Hide only its age notice
# until the base image is rebuilt with current browser data.
ENV BROWSERSLIST_IGNORE_OLD_DATA=true

# Copy app source code
COPY . .

#Expose port and start application
EXPOSE 3000
CMD ["npm", "start"]
