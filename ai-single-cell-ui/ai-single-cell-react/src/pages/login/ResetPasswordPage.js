import React, { useState } from 'react';
import { useParams } from 'react-router-dom';

const ResetPasswordPage = () => {
  const { token } = useParams();
  const [newPassword, setNewPassword] = useState('');
  const [confirmPassword, setConfirmPassword] = useState('');
  const [message, setMessage] = useState('');

  const handleSubmit = (e) => {
    console.log(token);
    e.preventDefault();
    if (newPassword !== confirmPassword) {
      setMessage('Passwords do not match');
      return;
    }

    fetch(`http://${process.env.REACT_APP_HOST_URL}:3001/api/reset-password`, {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({ token, newPassword,confirmPassword }),
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
      <h1>Reset Password</h1>
      <form onSubmit={handleSubmit}>
        <input
          type='password'
          value={newPassword}
          onChange={(e) => setNewPassword(e.target.value)}
          placeholder='Enter new password'
          required
        />
        <input
          type='password'
          value={confirmPassword}
          onChange={(e) => setConfirmPassword(e.target.value)}
          placeholder='Confirm new password'
          required
        />
        <button type='submit'>Reset Password</button>
      </form>
      {message && <p>{message}</p>}
    </div>
  );
};

export default ResetPasswordPage;
