import React, { useState } from 'react';

const TextWithEllipsis = ({ text, maxLength }) => {

  const displayText = text.length > maxLength ? text.slice(0, maxLength) : text;

  return (
    <div className="text-container">
        <div className="text-content" title={text}>
            {displayText}{text.length > maxLength && <div className="ellipsis">...</div>}

        </div>
    </div>
  );
};

export default TextWithEllipsis;
