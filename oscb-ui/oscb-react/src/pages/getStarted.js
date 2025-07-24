import { useState, useEffect } from "react"
import axios from 'axios';
import ReactMarkdown from "react-markdown"
import 'github-markdown-css';
// import LeftNav from "../components/LeftNavigation/leftNav";
import gfm from "remark-gfm";
// import rehypeRaw from 'rehype-raw'
import remarkImgToJsx from "remark-unwrap-images";
// import rehypeSanitize from 'rehype-sanitize'
import RightRail from "../components/RightNavigation/rightRail";
import { DIRECTUS_URL } from '../constants/declarations'
import { getCookie } from "../utils/utilFunctions";


export default function GetStarted() {
    const [markdownText,setMarkdownText] = useState('')
    let jwtToken = getCookie('jwtToken');

    useEffect(() => {
        async function fetchFileData() {
        try {
            const response = await axios.get(DIRECTUS_URL + "/items/filemappings?filter[filename]=getstarted");
            const data = response.data.data;
            
            if(data.length === 1) {
                const fileMappingObject = data[0];
                const fileID = fileMappingObject.fileID;
                if(fileID !== null) {
                    fetch(DIRECTUS_URL + "/assets/" + fileID)
                    .then(response => response.text())
                    .then(data => setMarkdownText(data))
                    .catch(error => console.error('Error retrieving markdown:', error));
                }
            }
            
          } catch (error) {
            console.error('Error retrieving data:', error);
          }
        }
        
        fetchFileData();
      }, []);
    
    return(
        
        <div className="getStarted-container">

            <div className="left-nav">
                {/* <LeftNav /> */}
            </div>

            <div className={(jwtToken === undefined || jwtToken === '') ? 'main-content-no-scroll' : 'main-content'}>
                <ReactMarkdown plugins={[gfm, remarkImgToJsx]} children={markdownText}/>
            </div>

            <div className="right-rail">
                <RightRail />
            </div>

            <div style={{display: 'none'}}>
                <script type='text/javascript' id='clustrmaps' src='//cdn.clustrmaps.com/map_v2.js?cl=fcfafa&w=a&t=tt&d=lOyy3dFp22wbqbeXPEE1e2nJSb_u_4KqYJPohHA8M4I&co=d8e1e8&cmo=f4acba&cmn=ff5353&ct=494545'></script>
            </div>  
            
        </div>
    )
}