import React from 'react';
import JobPlot from './JobPlot';

const Footer = () => {
    return (
        <footer style={{ borderTop: '1px solid #e5e7eb', backgroundColor: '#fff', marginTop: '24px' }}>
            {/* Job Plot Section */}
            <div style={{ maxWidth: '1536px', margin: '0 auto', padding: '12px 16px' }}>
                <JobPlot />
            </div>

            {/* Citation Section — centered */}
            <div style={{ maxWidth: '900px', margin: '0 auto', padding: '24px 16px 16px', textAlign: 'center' }}>
                <h3 style={{ fontSize: '16px', fontWeight: 700, color: '#1e3a5f', margin: '0 0 10px' }}>
                    Cite Single-Cell.AI
                </h3>
                <p style={{ fontSize: '13px', color: '#444', margin: '0 0 4px', lineHeight: 1.6 }}>
                    The Single-Cell.AI Community. "The platform for accessible, reproducible, and collaborative data analyses: 2024 update."
                </p>
                <p style={{ fontSize: '13px', color: '#666', fontStyle: 'italic', margin: 0, lineHeight: 1.6 }}>
                    Nucleic Acids Res. 2024, 52(W1):W83-W94. doi:10.1093/nar/gkae410 · <a href="#" style={{ color: '#2563eb', textDecoration: 'none' }}>full citation guide</a>
                </p>
            </div>

            {/* Horizontal Divider */}
            <hr style={{ border: 'none', borderTop: '1px solid #e0e0e0', margin: '0 auto', maxWidth: '900px' }} />

            {/* Institutional Logos Row — centered */}
            <div style={{
                maxWidth: '900px',
                margin: '0 auto',
                padding: '20px 16px',
                display: 'flex',
                flexWrap: 'wrap',
                justifyContent: 'center',
                alignItems: 'center',
                gap: '32px',
            }}>
                <div style={{ display: 'flex', alignItems: 'center', gap: '6px' }}>
                    <div style={{
                        width: '28px', height: '28px', backgroundColor: '#1e3a8a', color: '#fff',
                        borderRadius: '4px', display: 'flex', alignItems: 'center', justifyContent: 'center',
                        fontSize: '10px', fontWeight: 800,
                    }}>PS</div>
                    <span style={{ fontSize: '16px', fontWeight: 700, color: '#1e3a8a' }}>PennState</span>
                </div>
                <div style={{ display: 'flex', alignItems: 'center', gap: '6px' }}>
                    <div style={{
                        width: '28px', height: '28px', backgroundColor: '#1e3a5f', color: '#fff',
                        borderRadius: '4px', display: 'flex', alignItems: 'center', justifyContent: 'center',
                        fontSize: '8px', fontWeight: 800, fontFamily: 'serif',
                    }}>JHU</div>
                    <span style={{ fontSize: '16px', fontWeight: 700, color: '#1e3a5f', fontFamily: 'serif' }}>Johns Hopkins</span>
                </div>
                <span style={{ fontSize: '22px', fontWeight: 900, color: '#c41230', letterSpacing: '-0.5px' }}>TACC</span>
                <span style={{ fontSize: '22px', fontWeight: 900, color: '#0d9488', letterSpacing: '-0.5px' }}>ACCESS</span>
                <span style={{ fontSize: '22px', fontWeight: 700, color: '#b91c1c', fontStyle: 'italic' }}>Jetstream2</span>
            </div>

            {/* Two-Column Description — side by side */}
            <div style={{
                maxWidth: '900px',
                margin: '0 auto',
                padding: '0 16px 20px',
                display: 'grid',
                gridTemplateColumns: '1fr 1fr',
                gap: '24px',
            }}>
                <p style={{ fontSize: '12px', color: '#555', lineHeight: 1.6, margin: 0 }}>
                    The Single-Cell.AI Team is a part of the Center for Comparative Genomics and Bioinformatics at Penn State and the Department of Biology at Johns Hopkins University.
                </p>
                <p style={{ fontSize: '12px', color: '#555', lineHeight: 1.6, margin: 0 }}>
                    This instance of Single-Cell.AI is utilizing infrastructure generously provided by the Texas Advanced Computing Center. Additional resources are provided primarily on the Jetstream2 cloud via ACCESS, and with support from the National Science Foundation.
                </p>
            </div>

            {/* Bottom Version Line — centered */}
            <div style={{
                borderTop: '1px solid #eee',
                padding: '10px 16px',
                textAlign: 'center',
                fontSize: '11px',
                color: '#999',
            }}>
                Single-Cell.AI version 26.0.1.dev1, commit <span style={{ fontFamily: 'monospace', fontSize: '10px' }}>a4a38188871c9dad2b14b5a6a53f803a4a3a56e4</span>
            </div>
        </footer>
    );
};

export default Footer;
