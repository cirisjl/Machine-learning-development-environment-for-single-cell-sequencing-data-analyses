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
      "Disease Status (Specimen)": {
        "ui:placeholder": "Select an option",
        'ui:widget': () => (
          <div className='common-row-wrap'>
            <MyCreatableSelect fieldName="Disease Status (Specimen)" options={commonOptions} />
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
      "Sample Type": {
        "ui:placeholder": "Select an option",
        'ui:widget': () => (
          <div className='common-row-wrap'>
            <MyCreatableSelect fieldName="Sample Type" options={commonOptions}/>
          </div>
        ),
      },
      "Anatomical Entity": {
        "ui:placeholder": "Select an option",
        'ui:widget': () => (
          <div className='common-row-wrap'>
            <MyCreatableSelect fieldName="Anatomical Entity" options={commonOptions}/>
          </div>
        ),
      },
      "Model Organ": {
        "ui:placeholder": "Select an option",
        'ui:widget': () => (
          <div className='common-row-wrap'>
            <MyCreatableSelect fieldName="Model Organ" options={commonOptions}/>
          </div>
        ),
      },
      "Library Construction Method": {
        "ui:placeholder": "Select an option",
        'ui:widget': () => (
          <div className='common-row-wrap'>
            <MyCreatableSelect fieldName="Library Construction Method" options={commonOptions}/>
          </div>
        ),
      },
      "Nucleic Acid Source": {
        "ui:placeholder": "Select an option",
        'ui:widget': () => (
          <div className='common-row-wrap'>
            <MyCreatableSelect fieldName="Nucleic Acid Source" options={commonOptions}/>
          </div>
        ),
      },
      "Development Stage": {
        "ui:placeholder": "Select an option",
        'ui:widget': () => (
          <div className='common-row-wrap'>
            <MyCreatableSelect fieldName="Development Stage" options={commonOptions}/>
          </div>
        ),
      },
      "Source": {
        "ui:placeholder": "Select an option",
        'ui:widget': () => (
          <div className='common-row-wrap'>
            <MyCreatableSelect fieldName="Source" options={commonOptions}/>
          </div>
        ),
      },
};
};

export default createUISchema;
