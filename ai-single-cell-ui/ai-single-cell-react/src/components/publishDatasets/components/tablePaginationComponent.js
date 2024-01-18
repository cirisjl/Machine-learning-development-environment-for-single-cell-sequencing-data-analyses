import { FontAwesomeIcon } from '@fortawesome/react-fontawesome';
import React from 'react';

const Pagination = ({ pagination, onPageChange }) => {
    const { page, pageCount } = pagination;

    const goToPage = newPage => {
        if (newPage >= 1 && newPage <= pageCount) {
            onPageChange(newPage);
        }
    };

    return (
        <div className="pagination">
            <button onClick={() => goToPage(1)} disabled={page === 1}> <FontAwesomeIcon icon={'angles-left'} /></button>
            <button onClick={() => goToPage(page - 1)} disabled={page === 1}><FontAwesomeIcon icon={'angle-left'} /></button>
            <span>Page {page} of {pageCount}</span>
            <button onClick={() => goToPage(page + 1)} disabled={page === pageCount}><FontAwesomeIcon icon={'angle-right'} /></button>
            <button onClick={() => goToPage(pageCount)} disabled={page === pageCount}><FontAwesomeIcon icon={'angles-right'} /></button>
        </div>
    );
};

export default Pagination