import React , {useState, useEffect} from 'react';

function AlertMessageComponent({ message, setHasMessage, setMessage }) {

useEffect(() => {
    const timeoutId = setTimeout(() => {
        setMessage('');
        setHasMessage(false);
    }, 5000);
    // Return a cleanup function to cancel the timeout when the component unmounts
    return () => clearTimeout(timeoutId);
    }, [message]);


  return (
    <div>
        {message && 
        <div className='message-box' style={{ backgroundColor: '#bdf0c0' }}>
          <div style={{ textAlign: 'center' }}>
            <p>{message}</p>
          </div>
        </div>
        }       
    </div>
  );
}

export default AlertMessageComponent;
