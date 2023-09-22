//  // uiSchema.js

//  import MyCreatableSelect from "./Components/Creatable";

//  // Define the function that takes commonOptions as an argument
//  const createUISchema = (commonOptions) => {
//    return {
//      "parameters": {
//          "dataset": {
//            "ui:widget": "textarea",
//            "ui:placeholder": "Enter Text"
//          },
//          "submission_date": {
//            "classNames": "sub-category",
//          }
//      }
//    };
//  };
 
//  export default createUISchema;

export const uiSchema = {
  "parameters": {
      "dataset": {
        "ui:widget": "textarea",
        "ui:placeholder": "Enter Text"
      }
  }
};