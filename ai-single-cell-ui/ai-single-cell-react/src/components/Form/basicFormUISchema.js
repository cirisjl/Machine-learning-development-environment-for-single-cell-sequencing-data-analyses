  // uiSchema.js

import MyCreatableSelect from "./Components/Creatable";

// Define the function that takes commonOptions as an argument
const createUISchema = (commonOptions) => {
  return {
    "Dataset": {
      "ui:placeholder": "Enter the Dataset name"
    },
    "Submission Date": {
      "className":"date"
    },
};
};

export default createUISchema;
