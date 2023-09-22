import MyCreatableSelect from "../../components/Form/Components/Creatable";

  // Define the function that takes commonOptions as an argument
  const createUISchema = (commonOptions) => {
    return {
      "Dataset": {
        "ui:placeholder": "Enter the Dataset name"
      },
      "Task": {
        "ui:placeholder": "Select/Create an Option",
        'ui:widget': () => (
          <div className='common-row-wrap'>
            <MyCreatableSelect fieldName="Task" options={commonOptions}/>
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
              <MyCreatableSelect fieldName="Author" options={commonOptions} />
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
              <MyCreatableSelect fieldName="Species" options={commonOptions} />
            </div>
          ),
      },
      "Sample Type": {
        "ui:placeholder": "Select/Create an Option",
        'ui:widget': () => (
          <div className='common-row-wrap'>
            <MyCreatableSelect fieldName="Sample Type" options={commonOptions}/>
          </div>
        ),
      },
      "Anatomical Entity": {
        "ui:placeholder": "Select/Create an Option",
        'ui:widget': () => (
          <div className='common-row-wrap'>
            <MyCreatableSelect fieldName="Anatomical Entity" options={commonOptions}/>
          </div>
        ),
      },
      "Organ Part": {
        "ui:placeholder": "Select/Create an Option",
        'ui:widget': () => (
          <div className='common-row-wrap'>
            <MyCreatableSelect fieldName="Organ Part" options={commonOptions}/>
          </div>
        ),
      },
      "Model Organ": {
        "ui:placeholder": "Select/Create an Option",
        'ui:widget': () => (
          <div className='common-row-wrap'>
            <MyCreatableSelect fieldName="Model Organ" options={commonOptions}/>
          </div>
        ),
      },
      "Selected Cell Types": {
        "ui:placeholder": "Select/Create an Option",
        'ui:widget': () => (
          <div className='common-row-wrap'>
            <MyCreatableSelect fieldName="Selected Cell Types" options={commonOptions}/>
          </div>
        ),
      },
      
      "Library Construction Method": {
        "ui:placeholder": "Select/Create an Option",
        'ui:widget': () => (
          <div className='common-row-wrap'>
            <MyCreatableSelect fieldName="Library Construction Method" options={commonOptions}/>
          </div>
        ),
      },
      
      "Nucleic Acid Source": {
        "ui:placeholder": "Select/Create an Option",
        'ui:widget': () => (
          <div className='common-row-wrap'>
            <MyCreatableSelect fieldName="Nucleic Acid Source" options={commonOptions}/>
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
            <MyCreatableSelect fieldName="Disease Status (Specimen)" options={commonOptions} />
          </div>
        ),
      },
      
      "Disease Status (Donor)": {
        "ui:placeholder": "Select/Create an Option",
        'ui:widget': () => (
          <div className='common-row-wrap'>
            <MyCreatableSelect fieldName="Disease Status (Donor)" options={commonOptions} />
          </div>
        ),
      },
      
      "Development Stage": {
        "ui:placeholder": "Select/Create an Option",
        'ui:widget': () => (
          <div className='common-row-wrap'>
            <MyCreatableSelect fieldName="Development Stage" options={commonOptions}/>
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
            <MyCreatableSelect fieldName="Cell Count Estimate" options={commonOptions}/>
          </div>
        ),
      },
      "Source": {
        "ui:placeholder": "Select/Create an Option",
        'ui:widget': () => (
          <div className='common-row-wrap'>
            <MyCreatableSelect fieldName="Source" options={commonOptions}/>
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
  };
  
  export default createUISchema;
  