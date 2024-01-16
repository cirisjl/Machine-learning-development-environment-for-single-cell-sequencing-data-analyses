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
            <button onClick={() => goToPage(1)} disabled={page === 1}>First</button>
            <button onClick={() => goToPage(page - 1)} disabled={page === 1}>Previous</button>
            <span>Page {page} of {pageCount}</span>
            <button onClick={() => goToPage(page + 1)} disabled={page === pageCount}>Next</button>
            <button onClick={() => goToPage(pageCount)} disabled={page === pageCount}>Last</button>
        </div>
    );
};

export default Pagination