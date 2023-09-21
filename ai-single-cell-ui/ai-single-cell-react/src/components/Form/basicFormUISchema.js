// import MyCreatableSelect from "./Components/Creatable";
// const uiSchema = {
//         "Task": {
//           "ui:placeholder": "Select an option",
//           'ui:widget': () => (
//             <div className='common-row-wrap'>
//               <MyCreatableSelect fieldName="task"/>
//             </div>
//           ),
//         },
//         "Author": {
//             "ui:placeholder": "Select an option",
//             'ui:widget': () => (
//               <div className='common-row-wrap'>
//                 <MyCreatableSelect fieldName="author" />
//               </div>
//             ),
//         },
//         "Species": {
//             "ui:placeholder": "Select an option",
//             'ui:widget': () => (
//               <div className='common-row-wrap'>
//                 <MyCreatableSelect fieldName="species" />
//               </div>
//             ),
//         },
//         "Organ Part": {
//             "ui:placeholder": "Select an option",
//             'ui:widget': () => (
//               <div className='common-row-wrap'>
//                 <MyCreatableSelect fieldName="organ-part" />
//               </div>
//             ),
//           },
//           "Selected Cell Types": {
//             "ui:placeholder": "Select an option",
//             'ui:widget': () => (
//               <div className='common-row-wrap'>
//                 <MyCreatableSelect fieldName="selected-cell-types" />
//               </div>
//             ),
//           },
//           "Disease Status (Donor)": {
//             "ui:placeholder": "Select an option",
//             'ui:widget': () => (
//               <div className='common-row-wrap'>
//                 <MyCreatableSelect fieldName="disease-status-donor" />
//               </div>
//             ),
//           },
//           "Cell Count Estimate": {
//             "ui:placeholder": "Select an option",
//             'ui:widget': () => (
//               <div className='common-row-wrap'>
//                 <MyCreatableSelect fieldName="cell-count-estimate" />
//               </div>
//             ),
//           },
//   };
  

//   export default uiSchema;


  // uiSchema.js

import MyCreatableSelect from "./Components/Creatable";

// Define the function that takes commonOptions as an argument
const createUISchema = (commonOptions) => {
  return {
    "Task": {
      "ui:placeholder": "Select an option",
      'ui:widget': () => (
        <div className='common-row-wrap'>
          <MyCreatableSelect fieldName="Task" options={commonOptions}/>
        </div>
      ),
    },
    "Author": {
        "ui:placeholder": "Select an option",
        'ui:widget': () => (
          <div className='common-row-wrap'>
            <MyCreatableSelect fieldName="Author" options={commonOptions} />
          </div>
        ),
    },
    "Species": {
        "ui:placeholder": "Select an option",
        'ui:widget': () => (
          <div className='common-row-wrap'>
            <MyCreatableSelect fieldName="Species" options={commonOptions} />
          </div>
        ),
    },
    "Organ Part": {
        "ui:placeholder": "Select an option",
        'ui:widget': () => (
          <div className='common-row-wrap'>
            <MyCreatableSelect fieldName="Organ Part" options={commonOptions}/>
          </div>
        ),
      },
      "Selected Cell Types": {
        "ui:placeholder": "Select an option",
        'ui:widget': () => (
          <div className='common-row-wrap'>
            <MyCreatableSelect fieldName="Selected Cell Types" options={commonOptions}/>
          </div>
        ),
      },
      "Disease Status (Donor)": {
        "ui:placeholder": "Select an option",
        'ui:widget': () => (
          <div className='common-row-wrap'>
            <MyCreatableSelect fieldName="Disease Status (Donor)" options={commonOptions} />
          </div>
        ),
      },
      "Cell Count Estimate": {
        "ui:placeholder": "Select an option",
        'ui:widget': () => (
          <div className='common-row-wrap'>
            <MyCreatableSelect fieldName="Cell Count Estimate" options={commonOptions}/>
          </div>
        ),
      },
};
};

export default createUISchema;
