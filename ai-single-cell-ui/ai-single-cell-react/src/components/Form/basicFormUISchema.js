import MyCreatableSelect from "./Components/Creatable";
const uiSchema = {
        "Task": {
          "ui:placeholder": "Select an option",
          'ui:widget': () => (
            <div className='common-row-wrap'>
              <MyCreatableSelect fieldName="task"/>
            </div>
          ),
        },
        "Author": {
            "ui:placeholder": "Select an option",
            'ui:widget': () => (
              <div className='common-row-wrap'>
                <MyCreatableSelect fieldName="author" />
              </div>
            ),
        },
        "Species": {
            "ui:placeholder": "Select an option",
            'ui:widget': () => (
              <div className='common-row-wrap'>
                <MyCreatableSelect fieldName="species" />
              </div>
            ),
        },
        "Organ Part": {
            "ui:placeholder": "Select an option",
            'ui:widget': () => (
              <div className='common-row-wrap'>
                <MyCreatableSelect fieldName="organ-part" />
              </div>
            ),
          },
          "Selected Cell Types": {
            "ui:placeholder": "Select an option",
            'ui:widget': () => (
              <div className='common-row-wrap'>
                <MyCreatableSelect fieldName="selected-cell-types" />
              </div>
            ),
          },
          "Disease Status (Donor)": {
            "ui:placeholder": "Select an option",
            'ui:widget': () => (
              <div className='common-row-wrap'>
                <MyCreatableSelect fieldName="disease-status-donor" />
              </div>
            ),
          },
          "Cell Count Estimate": {
            "ui:placeholder": "Select an option",
            'ui:widget': () => (
              <div className='common-row-wrap'>
                <MyCreatableSelect fieldName="cell-count-estimate" />
              </div>
            ),
          },
  };
  

  export default uiSchema;