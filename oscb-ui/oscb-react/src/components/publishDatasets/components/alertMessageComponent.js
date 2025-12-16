import React , {useState, useEffect} from 'react';

function AlertMessageComponent({ message, setHasMessage, setMessage, isError }) {

  const backgroundColor = isError ? '#f0c0c0' : '#bdf0c0';

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
        <div className='message-box' style={{ backgroundColor: backgroundColor}}>
          <div style={{ textAlign: 'center' }}>
            <p>{message}</p>
          </div>
        </div>
        }       
    </div>
  );
}

export default AlertMessageComponent;
