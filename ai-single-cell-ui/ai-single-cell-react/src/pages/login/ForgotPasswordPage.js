import React, { useState } from 'react';

const ForgotPasswordPage = () => {
  const [email, setEmail] = useState('');
  const [message, setMessage] = useState('');

  const handleEmailChange = (e) => {
    setEmail(e.target.value);
  };

  const handleSubmit = (e) => {
    e.preventDefault();

    fetch(`http://${process.env.REACT_APP_HOST_URL}:3001/api/forgot-password`, {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({ email }),
    })
    .then(response => {
      if (!response.ok) {
        return response.json().then(error => { throw new Error(error.message) });
      }
      return response.json();
    })
    .then(data => {
      setMessage(data.message);
    })
    .catch(error => {
      console.error('Error occurred:', error);
      setMessage('An error occurred, please try again.');
    });
  };

  return (
    <div>
      <h1>Forgot Password</h1>
      <form onSubmit={handleSubmit}>
        <input
          type='email'
          value={email}
          onChange={handleEmailChange}
          placeholder='Enter your email'
          required
        />
        <button type='submit'>Reset Password</button>
      </form>
      {message && <p>{message}</p>}
    </div>
  );
};

export default ForgotPasswordPage;