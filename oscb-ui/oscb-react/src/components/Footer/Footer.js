import React from 'react';
import JobPlot from './JobPlot';

const Footer = () => {
    return (
        <footer className="bg-white border-t border-gray-200 mt-10">
            {/* Plot Section */}
            <div className="mx-auto px-4 py-8 border-b border-gray-100" style={{ maxWidth: '1536px' }}>
                <JobPlot />
            </div>

            {/* Citation Section */}
            <div className="mx-auto px-4 py-10 flex flex-col items-center text-center" style={{ maxWidth: '1536px' }}>
                <h2 className="text-2xl font-bold text-gray-800 mb-4">Cite Single-Cell.AI</h2>
                <div className="max-w-4xl text-gray-600 mb-2">
                    The Single-Cell.AI Community. "The platform for accessible, reproducible, and collaborative data analyses: 2024 update."
                </div>
                <div className="max-w-4xl text-sm italic text-gray-500 mb-8">
                    Nucleic Acids Res. 2024, 52(W1):W83-W94. doi:10.1093/nar/gkae410 • full citation guide
                </div>
                
                {/* Institutional Logos / Placeholders */}
                <div className="flex flex-wrap justify-center items-center gap-6 mt-6 pb-6 opacity-90 transition-all">
                    <div className="flex items-center gap-2">
                        <div className="w-10 h-10 bg-blue-800 text-white rounded-md flex items-center justify-center font-bold text-xl border border-blue-900">PS</div>
                        <div className="text-xl font-bold text-blue-800">PennState</div>
                    </div>
                    <div className="flex items-center gap-2 px-6">
                        <div className="w-10 h-10 bg-blue-900 text-white rounded-md flex items-center justify-center font-serif font-bold text-xl border border-blue-950">JHU</div>
                        <div className="text-xl font-serif font-bold text-blue-900">Johns Hopkins</div>
                    </div>
                    <div className="text-3xl font-black text-blue-900 tracking-tighter px-4">TACC</div>
                    <div className="text-3xl font-black text-teal-600 tracking-tighter px-4">ACCESS</div>
                    <div className="text-3xl font-bold text-red-700 italic px-4">Jetstream2</div>
                </div>

                <div className="grid grid-cols-1 md:grid-cols-2 gap-8 text-left text-sm text-gray-600 mt-6 border-t border-gray-100 pt-8">
                    <div>
                        The Single-Cell.AI Team is a part of the Center for Comparative Genomics and Bioinformatics at Penn State and the Department of Biology at Johns Hopkins University.
                    </div>
                    <div>
                        This instance of Single-Cell.AI is utilizing infrastructure generously provided by the Texas Advanced Computing Center. Additional resources are provided primarily on the Jetstream2 cloud via ACCESS, and with support from the National Science Foundation.
                    </div>
                </div>
            </div>

            {/* Acknowledgement and Links */}
            <div className="bg-gray-100 py-4 px-4 text-sm text-gray-600 border-t border-gray-300">
                <div className="mx-auto flex flex-col md:flex-row justify-between items-center gap-4" style={{ maxWidth: '1536px' }}>
                    <div className="text-center md:text-left text-xs leading-relaxed max-w-2xl">
                        Single-Cell.AI is maintained largely by the Freiburg Galaxy Team but also collectively by groups and individuals from across Europe. All of the member sites in this repository contribute to the European Single-Cell.AI Project. For acknowledgement, please refer to the About section. All content on this site is available under CC0-1.0 unless otherwise specified.
                        <div className="mt-2 text-gray-400">
                            Single-Cell.AI version 26.0.1.dev1, commit a4a38188...
                        </div>
                    </div>
                    <div className="flex gap-4 items-center flex-wrap justify-center">
                        <a href="https://github.com" className="hover:text-indigo-600 transition-colors flex items-center gap-2">
                             <svg className="w-4 h-4 fill-current" viewBox="0 0 24 24"><path d="M12 .297c-6.63 0-12 5.373-12 12 0 5.303 3.438 9.8 8.205 11.385.6.113.82-.258.82-.577 0-.285-.01-1.04-.015-2.04-3.338.724-4.042-1.61-4.042-1.61C4.422 18.07 3.633 17.7 3.633 17.7c-1.087-.744.084-.729.084-.729 1.205.084 1.838 1.236 1.838 1.236 1.07 1.835 2.809 1.305 3.495.998.108-.776.417-1.305.76-1.605-2.665-.3-5.466-1.332-5.466-5.93 0-1.31.465-2.38 1.235-3.22-.135-.303-.54-1.523.105-3.176 0 0 1.005-.322 3.3 1.23.96-.267 1.98-.399 3-.405 1.02.006 2.04.138 3 .405 2.28-1.552 3.285-1.23 3.285-1.23.645 1.653.24 2.873.12 3.176.765.84 1.23 1.91 1.23 3.22 0 4.61-2.805 5.625-5.475 5.92.42.36.81 1.096.81 2.22 0 1.606-.015 2.896-.015 3.286 0 .315.21.69.825.57C20.565 22.092 24 17.592 24 12.297c0-6.627-5.373-12-12-12"/></svg>
                             Edit this page
                        </a>
                        <a href="mailto:contact@single-cell.ai" className="hover:text-indigo-600 transition-colors flex items-center gap-2">
                             <svg className="w-4 h-4 fill-current" viewBox="0 0 24 24"><path d="M20 4H4c-1.1 0-2 .9-2 2v12c0 1.1.9 2 2 2h16c1.1 0 2-.9 2-2V6c0-1.1-.9-2-2-2zm0 4l-8 5-8-5V6l8 5 8-5v2z"/></svg>
                             Contact Us
                        </a>
                        <a href="#" className="hover:text-indigo-600 transition-colors flex items-center gap-2">
                             <svg className="w-4 h-4 fill-current" viewBox="0 0 24 24"><path d="M4 11A9 9 0 0 1 13 20h3A12 12 0 0 0 4 8v3zM4 4V1A19 19 0 0 1 23 20h-3A16 16 0 0 0 4 4zm2 14a2 2 0 1 1-4 0 2 2 0 0 1 4 0z"/></svg>
                             Subscribe
                        </a>
                    </div>
                </div>
            </div>
        </footer>
    );
};

export default Footer;
