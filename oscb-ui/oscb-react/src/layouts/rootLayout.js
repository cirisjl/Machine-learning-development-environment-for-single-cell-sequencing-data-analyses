import { Outlet, Link, NavLink } from "react-router-dom"
import Authentication from "../components/Authentication/AuthForm";
// import Chatbot from "../components/RightNavigation/Chatbot";
// import SearchBox from "../components/Header/searchBar"
import React, { useState, useEffect } from "react";
import { deleteCookie, getCookie, isUserAuth } from "../utils/utilFunctions";
import { useNavigate } from 'react-router-dom';
import Footer from "../components/Footer/Footer";


export default function RootLayout() {

    const navigate = useNavigate();

    const [isLoginReq, setIsLoginReq] = useState(false);
    const [isUserLoggedIn, setIsUserLoggedIn] = useState(false);
    const [username, setUsername] = useState('');
    const [isAdmin, setIsAdmin] = useState(false);
    const [mobileMenuOpen, setMobileMenuOpen] = useState(false);

    const handleAuth = (event) => {
        event.preventDefault();
        setIsLoginReq(!isLoginReq);
    };
    const [hoveredChildIndex, setHoveredChildIndex] = useState(null);


    const handleMouseOver = (event) => {
        setHoveredChildIndex(parseInt(event.currentTarget.dataset.index));
    };

    const handleMouseOut = () => {
        setHoveredChildIndex(null);
    };

    const handleLogoutClick = () => {
        if (deleteCookie('jwtToken'))
            setIsUserLoggedIn(false);
        navigate('/getStarted');
        window.location.reload();

    }


    useEffect(() => {
        const jwtToken = getCookie('jwtToken');
        if (jwtToken) {
            // If the token exists, verify authenticity
            isUserAuth(jwtToken).then((authData) => {
                setIsUserLoggedIn(true);
                setUsername(authData.username);
                setIsAdmin(authData.isAdmin);
            })
        }
    }, [isUserLoggedIn]);

    // Shared icon style for consistency
    const iconStyle = { width: '16px', height: '16px', marginRight: '4px', flexShrink: 0 };

    return (
        <>
        <div className="auth-form-container">
            <Authentication isLoginReq={isLoginReq} handleAuth={handleAuth} />
        </div>
        <div style={{ width: '100%' }}>
                <header style={{ borderBottom: '1px solid #e5e7eb' }}>
                    <div style={{ padding: '0 16px', display: 'flex', height: '56px', alignItems: 'center', justifyContent: 'space-between' }}>
                        {/* Logo */}
                        <a style={{ display: 'flex', alignItems: 'center', gap: '8px', textDecoration: 'none', flexShrink: 0 }} href="/">
                            <img src={require("../assets/logo.png")} alt="Single-Cell.AI" style={{ width: '36px', height: '36px', objectFit: 'contain' }} />
                            <span style={{ fontWeight: 700, fontSize: '15px', color: '#111', whiteSpace: 'nowrap', letterSpacing: '-0.01em' }}>SINGLE-CELL.AI</span>
                        </a>

                        {/* Desktop Nav */}
                        <nav aria-label="Main" style={{ marginLeft: 'auto', display: 'none' }} className="desktop-nav">
                            <ul style={{ display: 'flex', alignItems: 'center', gap: '6px', listStyle: 'none', margin: 0, padding: 0, fontSize: '13px' }}>
                                <li data-index="0" onMouseOver={handleMouseOver} onMouseOut={handleMouseOut}>
                                    <NavLink to="getStarted" className="group" style={{ display: 'flex', alignItems: 'center', padding: '4px 6px', borderRadius: '4px', whiteSpace: 'nowrap', color: '#374151', textDecoration: 'none' }}>
                                        <svg style={iconStyle} xmlns="http://www.w3.org/2000/svg" viewBox="0 0 24 24">
                                            <path d="M20.23 7.24L12 12L3.77 7.24a1.98 1.98 0 0 1 .7-.71L11 2.76c.62-.35 1.38-.35 2 0l6.53 3.77c.29.173.531.418.7.71z" opacity=".25" fill="currentColor"></path>
                                            <path d="M12 12v9.5a2.09 2.09 0 0 1-.91-.21L4.5 17.48a2.003 2.003 0 0 1-1-1.73v-7.5a2.06 2.06 0 0 1 .27-1.01L12 12z" opacity=".5" fill="currentColor"></path>
                                            <path d="M20.5 8.25v7.5a2.003 2.003 0 0 1-1 1.73l-6.62 3.82c-.275.13-.576.198-.88.2V12l8.23-4.76c.175.308.268.656.27 1.01z" fill="currentColor"></path>
                                        </svg>
                                        Get Started
                                    </NavLink>
                                </li>
                                <li data-index="1" onMouseOver={handleMouseOver} onMouseOut={handleMouseOut}>
                                    <NavLink className="group" style={{ display: 'flex', alignItems: 'center', padding: '4px 6px', borderRadius: '4px', whiteSpace: 'nowrap', color: '#374151', textDecoration: 'none' }}>
                                        <svg style={iconStyle} xmlns="http://www.w3.org/2000/svg" viewBox="0 0 24 24">
                                            <path d="M15.273 18.728A6.728 6.728 0 1 1 22 11.999V12a6.735 6.735 0 0 1-6.727 6.728z" opacity=".5" fill="currentColor"></path>
                                            <path d="M8.727 18.728A6.728 6.728 0 1 1 15.455 12a6.735 6.735 0 0 1-6.728 6.728z" fill="currentColor"></path>
                                        </svg>
                                        Analyses
                                    </NavLink>
                                    <div className={hoveredChildIndex === 1 ? "suboptions-container" : "suboptions-container hide"}>
                                        <div className="rounded-xl border-gray-100 border styles-for-dropdown">
                                            <ul className="ul-suboptions">
                                                {isAdmin && (<li><NavLink to="manageOptions">Manage Form Options</NavLink></li>)}
                                                <li><NavLink to="mydata/upload-data">Upload Data</NavLink></li>
                                                <li><NavLink to="mydata">My Datasets</NavLink></li>
                                                <li><NavLink to="projectAdminPanel">My Projects</NavLink></li>
                                                <li><NavLink to="myJobs">My Jobs</NavLink></li>
                                                <li><NavLink to="mydata/workflows">Workflows</NavLink></li>
                                                <li><NavLink to="mydata/tools">Tools</NavLink></li>
                                            </ul>
                                        </div>
                                    </div>
                                </li>
                                <li data-index="2">
                                    <NavLink to="datasets" className="group" style={{ display: 'flex', alignItems: 'center', padding: '4px 6px', borderRadius: '4px', whiteSpace: 'nowrap', color: '#374151', textDecoration: 'none' }}>
                                        <svg style={iconStyle} xmlns="http://www.w3.org/2000/svg" viewBox="0 0 25 25">
                                            <ellipse cx="12.5" cy="5" fill="currentColor" fillOpacity="0.25" rx="7.5" ry="2"></ellipse>
                                            <path d="M12.5 15C16.6421 15 20 14.1046 20 13V20C20 21.1046 16.6421 22 12.5 22C8.35786 22 5 21.1046 5 20V13C5 14.1046 8.35786 15 12.5 15Z" fill="currentColor" opacity="0.5"></path>
                                            <path d="M12.5 7C16.6421 7 20 6.10457 20 5V11.5C20 12.6046 16.6421 13.5 12.5 13.5C8.35786 13.5 5 12.6046 5 11.5V5C5 6.10457 8.35786 7 12.5 7Z" fill="currentColor" opacity="0.5"></path>
                                            <path d="M5.23628 12C5.08204 12.1598 5 12.8273 5 13C5 14.1046 8.35786 15 12.5 15C16.6421 15 20 14.1046 20 13C20 12.8273 19.918 12.1598 19.7637 12C18.9311 12.8626 15.9947 13.5 12.5 13.5C9.0053 13.5 6.06886 12.8626 5.23628 12Z" fill="currentColor"></path>
                                        </svg>
                                        Datasets
                                    </NavLink>
                                </li>
                                <li data-index="3" onMouseOver={handleMouseOver} onMouseOut={handleMouseOut}>
                                    <NavLink className="group" style={{ display: 'flex', alignItems: 'center', padding: '4px 6px', borderRadius: '4px', whiteSpace: 'nowrap', color: '#374151', textDecoration: 'none' }}>
                                        <svg style={iconStyle} xmlns="http://www.w3.org/2000/svg" viewBox="0 0 24 24">
                                            <path d="M6 23H2a1 1 0 0 1-1-1v-8a1 1 0 0 1 1-1h4a1 1 0 0 1 1 1v8a1 1 0 0 1-1 1z" opacity=".25" fill="currentColor"></path>
                                            <path d="M14 23h-4a1 1 0 0 1-1-1V2a1 1 0 0 1 1-1h4a1 1 0 0 1 1 1v20a1 1 0 0 1-1 1z" fill="currentColor"></path>
                                            <path d="M22 23h-4a1 1 0 0 1-1-1V10a1 1 0 0 1 1-1h4a1 1 0 0 1 1 1v12a1 1 0 0 1-1 1z" opacity=".5" fill="currentColor"></path>
                                        </svg>
                                        Benchmarks
                                    </NavLink>
                                    <div className={hoveredChildIndex === 3 ? "suboptions-container" : "suboptions-container hide"}>
                                        <div className="rounded-xl border-gray-100 border styles-for-dropdown">
                                            <ul className="ul-suboptions">
                                                { /* <li><NavLink to="benchmarks">Overview</NavLink></li> */}
                                                <li><Link reloadDocument to="benchmarks/clustering">Clustering</Link></li>
                                                <li><Link reloadDocument to="benchmarks/imputation">Imputation</Link></li>
                                                <li><Link reloadDocument to="benchmarks/batch-integration">Batch Integration</Link></li>
                                                <li><Link reloadDocument to="benchmarks/multimodal-data-integration">Multimodal Data Integration</Link></li>
                                                <li><Link reloadDocument to="benchmarks/trajectory">Trajectory</Link></li>
                                                <li><Link reloadDocument to="benchmarks/cell-cell-communication">Cell-Cell Communication</Link></li>
                                                <li><Link reloadDocument to="benchmarks/cell-type-annotation">Cell Type Annotation</Link></li>
                                                {isAdmin && (<li><NavLink to="benchmarks/uploads">Create New Benchmarks</NavLink></li>)}
                                            </ul>
                                        </div>
                                    </div>
                                </li>
                                <li data-index="4">
                                    <NavLink to="updates" className="group" style={{ display: 'flex', alignItems: 'center', padding: '4px 6px', borderRadius: '4px', whiteSpace: 'nowrap', color: '#374151', textDecoration: 'none' }}>
                                        <svg style={iconStyle} xmlns="http://www.w3.org/2000/svg" fill="currentColor" viewBox="0 0 24 24">
                                            <path fillRule="evenodd" d="M3.559 4.544c.355-.35.834-.544 1.33-.544H19.11c.496 0 .975.194 1.33.544.356.35.559.829.559 1.331v9.25c0 .502-.203.981-.559 1.331-.355.35-.834.544-1.33.544H15.5l-2.7 3.6a1 1 0 0 1-1.6 0L8.5 17H4.889c-.496 0-.975-.194-1.33-.544A1.868 1.868 0 0 1 3 15.125v-9.25c0-.502.203-.981.559-1.331ZM7.556 7.5a1 1 0 1 0 0 2h8a1 1 0 0 0 0-2h-8Zm0 3.5a1 1 0 1 0 0 2H12a1 1 0 1 0 0-2H7.556Z" clipRule="evenodd" />
                                        </svg>
                                        Updates
                                    </NavLink>
                                </li>
                                <li data-index="5" onMouseOver={handleMouseOver} onMouseOut={handleMouseOut}>
                                    <NavLink to="doc" className="group" style={{ display: 'flex', alignItems: 'center', padding: '4px 6px', borderRadius: '4px', whiteSpace: 'nowrap', color: '#374151', textDecoration: 'none' }}>
                                        <svg style={iconStyle} xmlns="http://www.w3.org/2000/svg" fill="currentColor" viewBox="0 0 24 24">
                                            <path fillRule="evenodd" d="M11 4.717c-2.286-.58-4.16-.756-7.045-.71A1.99 1.99 0 0 0 2 6v11c0 1.133.934 2.022 2.044 2.007 2.759-.038 4.5.16 6.956.791V4.717Zm2 15.081c2.456-.631 4.198-.829 6.956-.791A2.013 2.013 0 0 0 22 16.999V6a1.99 1.99 0 0 0-1.955-1.993c-2.885-.046-4.76.13-7.045.71v15.081Z" clipRule="evenodd" />
                                        </svg>
                                        Docs
                                    </NavLink>
                                </li>
                                <li data-index="6">
                                    <NavLink to="team" className="group" style={{ display: 'flex', alignItems: 'center', padding: '4px 6px', borderRadius: '4px', whiteSpace: 'nowrap', color: '#374151', textDecoration: 'none' }}>
                                        <svg style={iconStyle} xmlns="http://www.w3.org/2000/svg" fill="currentColor" viewBox="0 0 24 24">
                                            <path fillRule="evenodd" d="M12 6a3.5 3.5 0 1 0 0 7 3.5 3.5 0 0 0 0-7Zm-1.5 8a4 4 0 0 0-4 4 2 2 0 0 0 2 2h7a2 2 0 0 0 2-2 4 4 0 0 0-4-4h-3Zm6.82-3.096a5.51 5.51 0 0 0-2.797-6.293 3.5 3.5 0 1 1 2.796 6.292ZM19.5 18h.5a2 2 0 0 0 2-2 4 4 0 0 0-4-4h-1.1a5.503 5.503 0 0 1-.471.762A5.998 5.998 0 0 1 19.5 18ZM4 7.5a3.5 3.5 0 0 1 5.477-2.889 5.5 5.5 0 0 0-2.796 6.293A3.501 3.501 0 0 1 4 7.5ZM7.1 12H6a4 4 0 0 0-4 4 2 2 0 0 0 2 2h.5a5.998 5.998 0 0 1 3.071-5.238A5.505 5.505 0 0 1 7.1 12Z" clipRule="evenodd" />
                                        </svg>
                                        Teams
                                    </NavLink>
                                </li>
                                <li data-index="7" onMouseOver={handleMouseOver} onMouseOut={handleMouseOut}>
                                    <NavLink className="group" style={{ display: 'flex', alignItems: 'center', padding: '4px 6px', borderRadius: '4px', whiteSpace: 'nowrap', color: '#374151', textDecoration: 'none' }}>
                                        <svg style={iconStyle} xmlns="http://www.w3.org/2000/svg" fill="currentColor" viewBox="0 0 24 24">
                                            <path fillRule="evenodd" d="M12 20a7.966 7.966 0 0 1-5.002-1.756l.002.001v-.683c0-1.794 1.492-3.25 3.333-3.25h3.334c1.84 0 3.333 1.456 3.333 3.25v.683A7.966 7.966 0 0 1 12 20ZM2 12C2 6.477 6.477 2 12 2s10 4.477 10 10c0 5.5-4.44 9.963-9.932 10h-.138C6.438 21.962 2 17.5 2 12Zm10-5c-1.84 0-3.333 1.455-3.333 3.25S10.159 13.5 12 13.5c1.84 0 3.333-1.455 3.333-3.25S13.841 7 12 7Z" clipRule="evenodd" />
                                        </svg>
                                        {isUserLoggedIn ? (<span>Hi, <strong>{username}</strong>!</span>) : (<span>Login/Sign Up</span>)}
                                    </NavLink>
                                    <div className={hoveredChildIndex === 7 ? "suboptions-container" : "suboptions-container hide"}>
                                        <div className="rounded-xl border-gray-100 border styles-for-dropdown">
                                            <ul className="ul-suboptions">
                                                {!isUserLoggedIn && (<li><NavLink to="SignUp">Sign Up</NavLink></li>)}
                                                {!isUserLoggedIn && (<li><NavLink to="forgot-password">Forgot Password</NavLink></li>)}
                                                {isUserLoggedIn && (<li><NavLink to="reset/:token">Reset Password</NavLink></li>)}
                                                {isUserLoggedIn ? (
                                                    <li><span style={{ cursor: 'pointer' }} onClick={handleLogoutClick}>Log Out</span></li>
                                                ) : (
                                                    <li><NavLink to="login">Log In</NavLink></li>
                                                )}
                                            </ul>
                                        </div>
                                    </div>
                                </li>
                            </ul>
                        </nav>

                        {/* Mobile Hamburger Button */}
                        <button
                            onClick={() => setMobileMenuOpen(!mobileMenuOpen)}
                            className="mobile-menu-btn"
                            aria-label="Toggle menu"
                            style={{
                                display: 'none', background: 'none', border: 'none', cursor: 'pointer',
                                padding: '8px', marginLeft: 'auto', fontSize: '20px', color: '#374151',
                            }}
                        >
                            {mobileMenuOpen ? '✕' : '☰'}
                        </button>
                    </div>

                    {/* Mobile Nav Overlay */}
                    {mobileMenuOpen && (
                        <div className="mobile-nav-overlay" style={{
                            borderTop: '1px solid #e5e7eb', backgroundColor: '#fff', padding: '8px 0',
                        }}>
                            <ul style={{ listStyle: 'none', margin: 0, padding: 0, fontSize: '14px' }}>
                                {[
                                    { to: 'getStarted', label: 'Get Started' },
                                    { to: 'mydata/upload-data', label: 'Upload Data' },
                                    { to: 'mydata', label: 'My Datasets' },
                                    { to: 'datasets', label: 'Datasets' },
                                    { to: 'benchmarks/clustering', label: 'Benchmarks' },
                                    { to: 'updates', label: 'Updates' },
                                    { to: 'doc', label: 'Docs' },
                                    { to: 'team', label: 'Teams' },
                                ].map((item) => (
                                    <li key={item.to}>
                                        <NavLink
                                            to={item.to}
                                            onClick={() => setMobileMenuOpen(false)}
                                            style={{ display: 'block', padding: '10px 20px', color: '#374151', textDecoration: 'none', borderBottom: '1px solid #f3f4f6' }}
                                        >
                                            {item.label}
                                        </NavLink>
                                    </li>
                                ))}
                                <li>
                                    {isUserLoggedIn ? (
                                        <span onClick={() => { handleLogoutClick(); setMobileMenuOpen(false); }} style={{ display: 'block', padding: '10px 20px', color: '#374151', cursor: 'pointer' }}>Log Out</span>
                                    ) : (
                                        <NavLink to="login" onClick={() => setMobileMenuOpen(false)} style={{ display: 'block', padding: '10px 20px', color: '#374151', textDecoration: 'none' }}>Login/Sign Up</NavLink>
                                    )}
                                </li>
                            </ul>
                        </div>
                    )}
                </header>
        </div>
        <div className="container">
            <div className="main-container">
                <main>
                    <Outlet isUserLoggedIn={isUserLoggedIn} />
                </main>
            </div>
        </div>
        <div>
            <Footer />
        </div>
        </>
    )
}
