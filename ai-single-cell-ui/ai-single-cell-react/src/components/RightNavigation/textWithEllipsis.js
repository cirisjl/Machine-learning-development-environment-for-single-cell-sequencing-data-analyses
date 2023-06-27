import React, { useState } from 'react';

const TextWithEllipsis = ({ text, maxLength }) => {
  const [isExpanded, setIsExpanded] = useState(false);

  const toggleExpand = () => {
    setIsExpanded(!isExpanded);
  };

  const displayText = isExpanded ? text : text.slice(0, maxLength);

  return (
    <div className="text-container">
      <div className="text-content" onMouseEnter={toggleExpand} onMouseLeave={toggleExpand}>
        {displayText}
      </div>
      {!isExpanded && text.length > maxLength && <div className="ellipsis">...</div>}
    </div>
  );
};

export default TextWithEllipsis;
