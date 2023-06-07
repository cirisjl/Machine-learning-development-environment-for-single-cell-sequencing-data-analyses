import {Outlet, NavLink} from "react-router-dom"
import Authentication from "../components/Authentication/AuthForm";
import SearchBox from "../components/Header/searchBar"
import React, { useState, useEffect } from "react";
import { deleteCookie, getCookie } from "../utils/utilFunctions";
import { useNavigate } from 'react-router-dom';


export default function RootLayout() {

    const navigate = useNavigate();

    const [isLoginReq, setIsLoginReq] = useState(false);

    const [isUserLoggedIn, setIsUserLoggedIn] = useState(false);

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
        if(deleteCookie('jwtToken'))
        setIsUserLoggedIn(false);
        navigate('/getStarted');
        window.location.reload();  
        
    }

    const handleLoginClick = () => {
        setIsUserLoggedIn(true);
    }

    useEffect(() => {
        const userCookie = getCookie('jwtToken');

        if(userCookie!== '') {
            setIsUserLoggedIn(true);
        }
        
      }, [isUserLoggedIn]);

    // useEffect(() => {
    //     const handleCookieChange = () => {
    //       const newCookieValue = document.cookie
    //         .split('; ')
    //         .find(row => row.startsWith('jwtToken='))
    //         ?.split('=')[1];
    
    //       setIsUserLoggedIn(Boolean(newCookieValue));
    //     };
    
    //     // Listen for changes to the cookie value
    //     window.addEventListener('change', handleCookieChange);
    
    //     // Clean up the event listener when the component unmounts
    //     return () => {
    //       window.removeEventListener('change', handleCookieChange);
    //     };
    //   }, []);
  
    return(
        <div className="container">
            <div className="auth-form-container">
                <Authentication isLoginReq={isLoginReq} handleAuth={handleAuth}/>
            </div>
            <div className ="header-container">
                <header className="border-b border-gray-100">
                    <div className="px-4 flex h-16 items-center">
                        <div className="flex flex-1 items-center">
                            <a className="mr-5 flex flex-none items-center lg:mr-6" href="/"><span className="hidden whitespace-nowrap text-lg font-bold md:block">OSCB</span></a>
                            <div className="relative flex-1 lg:max-w-sm mr-2 sm:mr-4 lg:mr-6">
                                <SearchBox />
                            </div>
                        </div>

                    {/* <div className="header-right"> */}
                        <nav aria-label="Main" className="ml-auto hidden lg:block">
                            <ul className="flex items-center space-x-2">
                                <li data-index="0" onMouseOver={handleMouseOver} onMouseOut={handleMouseOut}>
                                    <NavLink to="getStarted" className="group flex items-center py-0.5 dark:hover:text-gray-400 hover:text-indigo-700 whitespace-nowrap" >
                                        <svg className="mr-1 text-gray-400 group-hover:text-indigo-500" xmlns="http://www.w3.org/2000/svg" xmlnsXlink="http://www.w3.org/1999/xlink" aria-hidden="true" focusable="false" role="img" width="1em" height="1em" preserveAspectRatio="xMidYMid meet" viewBox="0 0 24 24">
                                            <path className="uim-quaternary" d="M20.23 7.24L12 12L3.77 7.24a1.98 1.98 0 0 1 .7-.71L11 2.76c.62-.35 1.38-.35 2 0l6.53 3.77c.29.173.531.418.7.71z" opacity=".25" fill="currentColor"></path>
                                            <path className="uim-tertiary" d="M12 12v9.5a2.09 2.09 0 0 1-.91-.21L4.5 17.48a2.003 2.003 0 0 1-1-1.73v-7.5a2.06 2.06 0 0 1 .27-1.01L12 12z" opacity=".5" fill="currentColor"></path>
                                            <path className="uim-primary" d="M20.5 8.25v7.5a2.003 2.003 0 0 1-1 1.73l-6.62 3.82c-.275.13-.576.198-.88.2V12l8.23-4.76c.175.308.268.656.27 1.01z" fill="currentColor"></path>
                                        </svg>
                                        Get Started
                                    </NavLink>
                                    <div className={hoveredChildIndex === 0 ? "suboptions-container" : "suboptions-container hide"} >
                                        <div className="rounded-xl border-gray-100 border styles-for-dropdown">
                                            <ul className="ul-suboptions">
                                                <li className="hover-link"><NavLink to="getStarted/benchmarks">Benchmarks</NavLink></li>
                                                <li><NavLink to="getStarted/mydata">My Data</NavLink></li>
                                            </ul>
                                        </div>
                                    </div>
                                   </li>
                                 <li data-index="1">
                                    <NavLink to="updates" className="group flex items-center py-0.5 dark:hover:text-gray-400 hover:text-indigo-700">
                                        <svg className="mr-1 text-gray-400 group-hover:text-red-500" xmlns="http://www.w3.org/2000/svg" xmlnsXlink="http://www.w3.org/1999/xlink" aria-hidden="true" focusable="false" role="img" width="1em" height="1em" preserveAspectRatio="xMidYMid meet" viewBox="0 0 25 25">
                                            <ellipse cx="12.5" cy="5" fill="currentColor" fillOpacity="0.25" rx="7.5" ry="2"></ellipse>
                                            <path d="M12.5 15C16.6421 15 20 14.1046 20 13V20C20 21.1046 16.6421 22 12.5 22C8.35786 22 5 21.1046 5 20V13C5 14.1046 8.35786 15 12.5 15Z" fill="currentColor" opacity="0.5"></path>
                                            <path d="M12.5 7C16.6421 7 20 6.10457 20 5V11.5C20 12.6046 16.6421 13.5 12.5 13.5C8.35786 13.5 5 12.6046 5 11.5V5C5 6.10457 8.35786 7 12.5 7Z" fill="currentColor" opacity="0.5"></path>
                                            <path d="M5.23628 12C5.08204 12.1598 5 12.8273 5 13C5 14.1046 8.35786 15 12.5 15C16.6421 15 20 14.1046 20 13C20 12.8273 19.918 12.1598 19.7637 12C18.9311 12.8626 15.9947 13.5 12.5 13.5C9.0053 13.5 6.06886 12.8626 5.23628 12Z" fill="currentColor"></path>
                                        </svg>
                                     Updates
                                    </NavLink>    
                                </li>
                                 <li data-index="2">
                                    <NavLink to="competitions" className="group flex items-center py-0.5 dark:hover:text-gray-400 hover:text-indigo-700">
                                        <svg className="mr-1 text-gray-400 group-hover:text-blue-500" xmlns="http://www.w3.org/2000/svg" xmlnsXlink="http://www.w3.org/1999/xlink" aria-hidden="true" focusable="false" role="img" width="1em" height="1em" viewBox="0 0 25 25">
                                            <path opacity=".5" d="M6.016 14.674v4.31h4.31v-4.31h-4.31ZM14.674 14.674v4.31h4.31v-4.31h-4.31ZM6.016 6.016v4.31h4.31v-4.31h-4.31Z" fill="currentColor"></path>
                                            <path opacity=".75" fillRule="evenodd" clipRule="evenodd" d="M3 4.914C3 3.857 3.857 3 4.914 3h6.514c.884 0 1.628.6 1.848 1.414a5.171 5.171 0 0 1 7.31 7.31c.815.22 1.414.964 1.414 1.848v6.514A1.914 1.914 0 0 1 20.086 22H4.914A1.914 1.914 0 0 1 3 20.086V4.914Zm3.016 1.102v4.31h4.31v-4.31h-4.31Zm0 12.968v-4.31h4.31v4.31h-4.31Zm8.658 0v-4.31h4.31v4.31h-4.31Zm0-10.813a2.155 2.155 0 1 1 4.31 0 2.155 2.155 0 0 1-4.31 0Z" fill="currentColor"></path>
                                            <path opacity=".25" d="M16.829 6.016a2.155 2.155 0 1 0 0 4.31 2.155 2.155 0 0 0 0-4.31Z" fill="currentColor"></path>
                                        </svg>
                                        Competitions
                                    </NavLink>
                                 </li>
                                 <li data-index="3" onMouseOver={handleMouseOver} onMouseOut={handleMouseOut}>
                                    <NavLink to="benchmarks" className="group flex items-center py-0.5 dark:hover:text-gray-400 hover:text-indigo-700">
                                        <svg xmlns="http://www.w3.org/2000/svg" xmlnsXlink="http://www.w3.org/1999/xlink" aria-hidden="true" role="img" className="mr-1 text-gray-400 group-hover:text-yellow-500" width="1em" height="1em" preserveAspectRatio="xMidYMid meet" viewBox="0 0 32 32">
                                            <path opacity="0.5" d="M20.9022 5.10334L10.8012 10.8791L7.76318 9.11193C8.07741 8.56791 8.5256 8.11332 9.06512 7.7914L15.9336 3.73907C17.0868 3.08811 18.5002 3.26422 19.6534 3.91519L19.3859 3.73911C19.9253 4.06087 20.5879 4.56025 20.9022 5.10334Z" fill="currentColor"></path>
                                            <path d="M10.7999 10.8792V28.5483C10.2136 28.5475 9.63494 28.4139 9.10745 28.1578C8.5429 27.8312 8.074 27.3621 7.74761 26.7975C7.42122 26.2327 7.24878 25.5923 7.24756 24.9402V10.9908C7.25062 10.3319 7.42358 9.68487 7.74973 9.1123L10.7999 10.8792Z" fill="currentColor" fillOpacity="0.75"></path>
                                            <path fillRule="evenodd" clipRule="evenodd" d="M21.3368 10.8499V6.918C21.3331 6.25959 21.16 5.61234 20.8346 5.03949L10.7971 10.8727L10.8046 10.874L21.3368 10.8499Z" fill="currentColor"></path>
                                            <path opacity="0.5" d="M21.7937 10.8488L10.7825 10.8741V28.5486L21.7937 28.5234C23.3344 28.5234 24.5835 27.2743 24.5835 25.7335V13.6387C24.5835 12.0979 23.4365 11.1233 21.7937 10.8488Z" fill="currentColor"></path>
                                        </svg>
                                        Benchmarks
                                    </NavLink>
                                    <div className={hoveredChildIndex === 3 ? "suboptions-container" : "suboptions-container hide"}>
                                        <div className="rounded-xl border-gray-100 border styles-for-dropdown">

                                        <ul className="ul-suboptions">
                                            <li><NavLink to="benchmarks/overview">Overview</NavLink></li>
                                            <li><NavLink to="benchmarks/clustering">Clustering</NavLink></li>
                                            <li><NavLink to="benchmarks/imputation">Imputation</NavLink></li>
                                            <li><NavLink to="benchmarks/maker-gene-identification">Marker Gene Identification</NavLink></li>
                                            <li><NavLink to="benchmarks/trajectory">Trajectory</NavLink></li>
                                            <li><NavLink to="benchmarks/cell-cell-communication">Cell-Cell Communication</NavLink></li>
                                            <li><NavLink to="benchmarks/multiomics-data-integration">Multiomics Data Integration</NavLink></li>
                                            <li><NavLink to="benchmarks/gene-regulatory-relations">Gene Regulatory Relations</NavLink></li>
                                            <li><NavLink to="benchmarks/genes-over-time">Genes Over Time</NavLink></li>
                                            <li><NavLink to="benchmarks/genes-over-condition">Genes Over Condition</NavLink></li>
                                            <li><NavLink to="benchmarks/cell-type">Cell Type</NavLink></li>
                                        </ul>
                                    </div>
                                    </div>
                                 </li>
                                 <li data-index="4" onMouseOver={handleMouseOver} onMouseOut={handleMouseOut}>
                                    <NavLink to="leaderboards" className="group flex items-center py-0.5 dark:hover:text-gray-400 hover:text-indigo-700">
                                        <svg className="text-gray-400 mr-1" xmlns="http://www.w3.org/2000/svg" xmlnsXlink="http://www.w3.org/1999/xlink" aria-hidden="true" focusable="false" role="img" width="1em" height="1em" preserveAspectRatio="xMidYMid meet" viewBox="0 0 24 24">
                                            <path className="uim-quaternary" d="M6 23H2a1 1 0 0 1-1-1v-8a1 1 0 0 1 1-1h4a1 1 0 0 1 1 1v8a1 1 0 0 1-1 1z" opacity=".25" fill="currentColor"></path>
                                            <path className="uim-primary" d="M14 23h-4a1 1 0 0 1-1-1V2a1 1 0 0 1 1-1h4a1 1 0 0 1 1 1v20a1 1 0 0 1-1 1z" fill="currentColor"></path>
                                            <path className="uim-tertiary" d="M22 23h-4a1 1 0 0 1-1-1V10a1 1 0 0 1 1-1h4a1 1 0 0 1 1 1v12a1 1 0 0 1-1 1z" opacity=".5" fill="currentColor"></path>
                                        </svg>
                                        Leaderboards
                                    </NavLink>
                                 </li>
                                 <li data-index="5" onMouseOver={handleMouseOver} onMouseOut={handleMouseOut}>
                                    <NavLink to="mydata" className="group flex items-center py-0.5 dark:hover:text-gray-400 hover:text-indigo-700">
                                        <svg className="text-gray-400 mr-1" xmlns="http://www.w3.org/2000/svg" xmlnsXlink="http://www.w3.org/1999/xlink" aria-hidden="true" role="img" width="1em" height="1em" preserveAspectRatio="xMidYMid meet" viewBox="0 0 24 24">
                                            <path className="uim-tertiary" d="M15.273 18.728A6.728 6.728 0 1 1 22 11.999V12a6.735 6.735 0 0 1-6.727 6.728z" opacity=".5" fill="currentColor"></path>
                                            <path className="uim-primary" d="M8.727 18.728A6.728 6.728 0 1 1 15.455 12a6.735 6.735 0 0 1-6.728 6.728z" fill="currentColor"></path>
                                        </svg>
                                        MyData
                                    </NavLink>
                                    <div className={hoveredChildIndex === 5 ? "suboptions-container" : "suboptions-container hide"}>
                                    <div className="rounded-xl border-gray-100 border styles-for-dropdown">
                                        <ul className="ul-suboptions">
                                            <li><NavLink to="mydata/upload-data">Upload Data</NavLink></li>
                                            <li><NavLink to="mydata/preview-datasets">My Datasets</NavLink></li>
                                            <li><NavLink to="mydata/workflows">Workflows</NavLink></li>
                                            <li><NavLink to="mydata/tools">Tools</NavLink></li>
                                        </ul>
                                    </div>
                                    </div>
                                 </li>
                                 <li data-index="6" onMouseOver={handleMouseOver} onMouseOut={handleMouseOut}>
                                    <NavLink to="docs" className="group flex items-center py-0.5 dark:hover:text-gray-400 hover:text-indigo-700">
                                        <svg className="mr-1 text-gray-400 group-hover:text-red-500" xmlns="http://www.w3.org/2000/svg" xmlnsXlink="http://www.w3.org/1999/xlink" aria-hidden="true" focusable="false" role="img" width="1em" height="1em" preserveAspectRatio="xMidYMid meet" viewBox="0 0 25 25">
                                            <ellipse cx="12.5" cy="5" fill="currentColor" fillOpacity="0.25" rx="7.5" ry="2"></ellipse>
                                            <path d="M12.5 15C16.6421 15 20 14.1046 20 13V20C20 21.1046 16.6421 22 12.5 22C8.35786 22 5 21.1046 5 20V13C5 14.1046 8.35786 15 12.5 15Z" fill="currentColor" opacity="0.5"></path>
                                            <path d="M12.5 7C16.6421 7 20 6.10457 20 5V11.5C20 12.6046 16.6421 13.5 12.5 13.5C8.35786 13.5 5 12.6046 5 11.5V5C5 6.10457 8.35786 7 12.5 7Z" fill="currentColor" opacity="0.5"></path>
                                            <path d="M5.23628 12C5.08204 12.1598 5 12.8273 5 13C5 14.1046 8.35786 15 12.5 15C16.6421 15 20 14.1046 20 13C20 12.8273 19.918 12.1598 19.7637 12C18.9311 12.8626 15.9947 13.5 12.5 13.5C9.0053 13.5 6.06886 12.8626 5.23628 12Z" fill="currentColor"></path>
                                        </svg>
                                        Docs
                                    </NavLink>
                                    <div className={hoveredChildIndex === 6 ? "suboptions-container" : "suboptions-container hide"}>
                                    <div className="rounded-xl border-gray-100 border styles-for-dropdown">
                                        <ul className="ul-suboptions">
                                            <li><NavLink to="docs/paper">Paper</NavLink></li>
                                            <li><NavLink to="docs/tutorial">Tutorial</NavLink></li>
                                        </ul>
                                    </div>
                                    </div>
                                 </li>
                                 <li data-index="7">
                                    <NavLink to="team" className="group flex items-center py-0.5 dark:hover:text-gray-400 hover:text-indigo-700">
                                        <svg xmlns="http://www.w3.org/2000/svg" xmlnsXlink="http://www.w3.org/1999/xlink" aria-hidden="true" role="img" className="mr-1 text-gray-400 group-hover:text-yellow-500" width="1em" height="1em" preserveAspectRatio="xMidYMid meet" viewBox="0 0 32 32">
                                            <path opacity="0.5" d="M20.9022 5.10334L10.8012 10.8791L7.76318 9.11193C8.07741 8.56791 8.5256 8.11332 9.06512 7.7914L15.9336 3.73907C17.0868 3.08811 18.5002 3.26422 19.6534 3.91519L19.3859 3.73911C19.9253 4.06087 20.5879 4.56025 20.9022 5.10334Z" fill="currentColor"></path>
                                            <path d="M10.7999 10.8792V28.5483C10.2136 28.5475 9.63494 28.4139 9.10745 28.1578C8.5429 27.8312 8.074 27.3621 7.74761 26.7975C7.42122 26.2327 7.24878 25.5923 7.24756 24.9402V10.9908C7.25062 10.3319 7.42358 9.68487 7.74973 9.1123L10.7999 10.8792Z" fill="currentColor" fillOpacity="0.75"></path>
                                            <path fillRule="evenodd" clipRule="evenodd" d="M21.3368 10.8499V6.918C21.3331 6.25959 21.16 5.61234 20.8346 5.03949L10.7971 10.8727L10.8046 10.874L21.3368 10.8499Z" fill="currentColor"></path>
                                            <path opacity="0.5" d="M21.7937 10.8488L10.7825 10.8741V28.5486L21.7937 28.5234C23.3344 28.5234 24.5835 27.2743 24.5835 25.7335V13.6387C24.5835 12.0979 23.4365 11.1233 21.7937 10.8488Z" fill="currentColor"></path>
                                        </svg>
                                        Teams
                                    </NavLink>
                                 </li>
                                 <li data-index="8" onMouseOver={handleMouseOver} onMouseOut={handleMouseOut}>
                                    <NavLink className="group flex items-center py-0.5 dark:hover:text-gray-400 hover:text-indigo-700">
                                        <svg className="mr-1 text-gray-400 group-hover:text-red-500" xmlns="http://www.w3.org/2000/svg" xmlnsXlink="http://www.w3.org/1999/xlink" aria-hidden="true" focusable="false" role="img" width="1em" height="1em" preserveAspectRatio="xMidYMid meet" viewBox="0 0 25 25">
                                            <ellipse cx="12.5" cy="5" fill="currentColor" fillOpacity="0.25" rx="7.5" ry="2"></ellipse>
                                            <path d="M12.5 15C16.6421 15 20 14.1046 20 13V20C20 21.1046 16.6421 22 12.5 22C8.35786 22 5 21.1046 5 20V13C5 14.1046 8.35786 15 12.5 15Z" fill="currentColor" opacity="0.5"></path>
                                            <path d="M12.5 7C16.6421 7 20 6.10457 20 5V11.5C20 12.6046 16.6421 13.5 12.5 13.5C8.35786 13.5 5 12.6046 5 11.5V5C5 6.10457 8.35786 7 12.5 7Z" fill="currentColor" opacity="0.5"></path>
                                            <path d="M5.23628 12C5.08204 12.1598 5 12.8273 5 13C5 14.1046 8.35786 15 12.5 15C16.6421 15 20 14.1046 20 13C20 12.8273 19.918 12.1598 19.7637 12C18.9311 12.8626 15.9947 13.5 12.5 13.5C9.0053 13.5 6.06886 12.8626 5.23628 12Z" fill="currentColor"></path>
                                        </svg>
                                        Login/SignUp
                                    </NavLink>
                                    <div className={hoveredChildIndex === 8 ? "suboptions-container" : "suboptions-container hide"}>
                                    <div className="rounded-xl border-gray-100 border styles-for-dropdown">
                                        <ul className="ul-suboptions">
                                            <li><NavLink to="login/my-projects">My Projects</NavLink></li>
                                            <li><NavLink to="login/my-data">My Data</NavLink></li>
                                            <li><NavLink to="login/reports">Reports</NavLink></li>
                                            <li><NavLink to="login/security">Security</NavLink></li>
                                            {isUserLoggedIn ? (
                                                <li><span style={{ cursor: 'pointer'}} onClick={handleLogoutClick}>Log Out</span></li>
                                            ) : (
                                                <li><NavLink to="login">Log In</NavLink></li>
                                            )}

                                        </ul>
                                    </div>
                                    </div>
                                 </li>
                             </ul>
                         </nav>
                     </div>
                 </header>
             </div>
             <div className="main-container">
                 <main>
                     <Outlet isUserLoggedIn={isUserLoggedIn}/>
                 </main>
             </div>
         </div>

    )
}
