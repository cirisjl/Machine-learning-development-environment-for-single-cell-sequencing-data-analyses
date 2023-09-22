  // uiSchema.js

  import MyCreatableSelect from "./Components/Creatable";

  // Define the function that takes commonOptions as an argument
  
  export const uiSchema = {
      "Dataset": {
        "ui:placeholder": "Enter the Dataset name"
      },
      "Task": {
        "ui:placeholder": "Select/Create an Option",
        'ui:widget': () => (
          <div className='common-row-wrap'>
            <MyCreatableSelect fieldName="Task" />
          </div>
        ),
      },
      "Downloads": {
        "ui:placeholder": "http://"
      },
      "Title": {
        "ui:placeholder": "Enter title"
      },
      "Author": {
          "ui:placeholder": "Select/Create an Option",
          'ui:widget': () => (
            <div className='common-row-wrap'>
              <MyCreatableSelect fieldName="Author" />
            </div>
          ),
      },
      "Reference (paper)": {
        "ui:placeholder": "Enter Reference"
      },
      "Abstract": {
        "ui:widget": "textarea"
      },
      "DOI": {
        "ui:widget": "textarea"
      },
      "Species": {
          "ui:placeholder": "Select/Create an Option",
          'ui:widget': () => (
            <div className='common-row-wrap'>
              <MyCreatableSelect fieldName="Species" />
            </div>
          ),
      },
      "Sample Type": {
        "ui:placeholder": "Select/Create an Option",
        'ui:widget': () => (
          <div className='common-row-wrap'>
            <MyCreatableSelect fieldName="Sample Type" />
          </div>
        ),
      },
      "Anatomical Entity": {
        "ui:placeholder": "Select/Create an Option",
        'ui:widget': () => (
          <div className='common-row-wrap'>
            <MyCreatableSelect fieldName="Anatomical Entity"/>
          </div>
        ),
      },
      "Organ Part": {
        "ui:placeholder": "Select/Create an Option",
        'ui:widget': () => (
          <div className='common-row-wrap'>
            <MyCreatableSelect fieldName="Organ Part"/>
          </div>
        ),
      },
      "Model Organ": {
        "ui:placeholder": "Select/Create an Option",
        'ui:widget': () => (
          <div className='common-row-wrap'>
            <MyCreatableSelect fieldName="Model Organ"/>
          </div>
        ),
      },
      "Selected Cell Types": {
        "ui:placeholder": "Select/Create an Option",
        'ui:widget': () => (
          <div className='common-row-wrap'>
            <MyCreatableSelect fieldName="Selected Cell Types"/>
          </div>
        ),
      },
      
      "Library Construction Method": {
        "ui:placeholder": "Select/Create an Option",
        'ui:widget': () => (
          <div className='common-row-wrap'>
            <MyCreatableSelect fieldName="Library Construction Method"/>
          </div>
        ),
      },
      
      "Nucleic Acid Source": {
        "ui:placeholder": "Select/Create an Option",
        'ui:widget': () => (
          <div className='common-row-wrap'>
            <MyCreatableSelect fieldName="Nucleic Acid Source"/>
          </div>
        ),
      },
  
      
      "Paired End": {
        "ui:widget": "radio" // could also be "select"
      },
      
      "Analysis Protocol": {
        "ui:placeholder": "Enter ...",
      },
      "Disease Status (Specimen)": {
        "ui:placeholder": "Select/Create an Option",
        'ui:widget': () => (
          <div className='common-row-wrap'>
            <MyCreatableSelect fieldName="Disease Status (Specimen)" />
          </div>
        ),
      },
      
      "Disease Status (Donor)": {
        "ui:placeholder": "Select/Create an Option",
        'ui:widget': () => (
          <div className='common-row-wrap'>
            <MyCreatableSelect fieldName="Disease Status (Donor)"/>
          </div>
        ),
      },
      
      "Development Stage": {
        "ui:placeholder": "Select/Create an Option",
        'ui:widget': () => (
          <div className='common-row-wrap'>
            <MyCreatableSelect fieldName="Development Stage"/>
          </div>
        ),
      },
      "Donor Count": {
        "ui:placeholder": "Enter ...",
      },
      "Cell Count Estimate": {
        "ui:placeholder": "Select/Create an Option",
        'ui:widget': () => (
          <div className='common-row-wrap'>
            <MyCreatableSelect fieldName="Cell Count Estimate"/>
          </div>
        ),
      },
      "Source": {
        "ui:placeholder": "Select/Create an Option",
        'ui:widget': () => (
          <div className='common-row-wrap'>
            <MyCreatableSelect fieldName="Source" />
          </div>
        ),
      },
      "Source Key": {
        "ui:placeholder": "Enter ...",
      },
      "Submission Date": {
        "className":"date"
      },
  };
  