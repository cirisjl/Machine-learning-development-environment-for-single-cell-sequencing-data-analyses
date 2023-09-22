 // uiSchema.js

 import MyCreatableSelect from "./Components/Creatable";

 // Define the function that takes commonOptions as an argument
 const createUISchema = (commonOptions) => {
   return {
     "parameters": {
         "dataset": {
           "ui:widget": "text-area",
           "ui:placeholder": "Enter Text"
         },
         "submission_date": {
           "classNames": "sub-category",
         }
     }
   };
 };
 
 export default createUISchema;