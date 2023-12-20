import React, { useState } from 'react';
import PublishDataset from '../publishDataset';

const FlowControl = () => {

  const [flow, setFlow] = useState('');


  const startFromBeginning = () => {
    setFlow('upload');
  };

  const startFromFourthStep = () => {
    setFlow('taskBuilder');
  };

  return (
    <div>
        {flow === '' && 
            <div className='flow-messaging'>
                <h3>If you want to create a new dataset and then build tasks, click on "Create a New Dataset". If you want to use existing datasets to build tasks, click on "Jump to Task Builder".</h3>
                <button onClick={startFromBeginning}>Create a new Dataset</button>
                <button onClick={startFromFourthStep}>Jump to Task builder</button>
            </div>
        }
      
        <div>
            {flow === 'upload' && <div><PublishDataset /></div>}
            {flow === 'taskBuilder' && <div>Step 4 Component or Content</div>}
        </div>
    </div>
  );
};

export default FlowControl;
