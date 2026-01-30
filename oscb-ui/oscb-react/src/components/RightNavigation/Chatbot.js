import React, { useState, useRef, useEffect } from 'react';
import axios from 'axios';
import { FontAwesomeIcon } from '@fortawesome/react-fontawesome';
import { faPaperPlane, faTrash, faUser, faRobot, faDna, faRotateRight, faMagic, faMinus, faExpand, faAnchor } from '@fortawesome/free-solid-svg-icons';
import { NODE_API_URL } from '../../constants/declarations';

import styled from 'styled-components';
import SingleCellLogo from '../../assets/single-cell-logo.png';

import ReactMarkdown from "react-markdown"
import 'github-markdown-css';
import { Prism as SyntaxHighlighter } from 'react-syntax-highlighter'
import { dark } from 'react-syntax-highlighter/dist/esm/styles/prism'
import remarkGfm from 'remark-gfm';
import rehypeRaw from 'rehype-raw'
import rehypeGithubAlerts from 'rehype-github-alert'
import { CopyToClipboard } from 'react-copy-to-clipboard';

// --- Styled Components ---

const Container = styled.div`
  display: flex;
  flex-direction: column;
  background-color: #ffffff;
  font-family: 'Inter', sans-serif;
  box-shadow: 0 4px 20px rgba(0, 0, 0, 0.15);
  position: fixed;
  bottom: 24px;
  right: 24px;
  z-index: 9999;
  border-radius: 16px;
  border: 1px solid #e2e8f0;
  overflow: hidden;
`;

const ResizeHandle = styled.div`
  position: absolute;
  top: 0;
  left: 0;
  width: 20px;
  height: 20px;
  cursor: nw-resize;
  z-index: 20;
  
  &::after {
    content: '';
    position: absolute;
    top: 6px;
    left: 6px;
    width: 6px;
    height: 6px;
    border-top: 2px solid #cbd5e1;
    border-left: 2px solid #cbd5e1;
  }

  &:hover::after {
    border-color: #64748b;
  }
`;

const Header = styled.div`
  padding: 16px 20px;
  border-bottom: 1px solid #f1f5f9;
  display: flex;
  justify-content: space-between;
  align-items: center;
  background-color: #ffffff;
  background-color: #ffffff;
  z-index: 10;
  cursor: grab;
  
  &:active {
    cursor: grabbing;
  }
`;

const MinimizedButton = styled.button`
  position: fixed;
  bottom: 24px;
  right: 24px;
  width: 60px;
  height: 60px;
  border-radius: 50%;
  background: transparent;
  color: white;
  border: none;
  /* box-shadow: 0 4px 12px rgba(13, 148, 136, 0.4); */
  cursor: pointer;
  z-index: 9999;
  display: flex;
  align-items: center;
  justify-content: center;
  transition: all 0.2s cubic-bezier(0.175, 0.885, 0.32, 1.275);

  &:hover {
    transform: scale(1.05);
  }
`;

const HeaderActions = styled.div`
  display: flex;
  align-items: center;
  gap: 8px;
`;

const HeaderTitleGroup = styled.div`
  display: flex;
  align-items: center;
  gap: 12px;
`;

const AvatarCircle = styled.div`
  width: 32px;
  height: 32px;
  border-radius: 50%;
  background: linear-gradient(135deg, #2dd4bf 0%, #0d9488 100%);
  display: flex;
  align-items: center;
  justify-content: center;
  color: white;
  box-shadow: 0 2px 4px rgba(13, 148, 136, 0.2);
`;

const TitleText = styled.h3`
  font-size: 14px;
  font-weight: 700;
  color: #1e293b;
  margin: 0;
  line-height: 1.2;
`;

const StatusIndicator = styled.div`
  display: flex;
  align-items: center;
  gap: 6px;
  margin-top: 2px;
`;

const StatusDot = styled.div`
  width: 6px;
  height: 6px;
  border-radius: 50%;
  background-color: #22c55e;
`;

const StatusText = styled.span`
  font-size: 10px;
  color: #94a3b8;
  font-weight: 500;
  text-transform: uppercase;
  letter-spacing: 0.5px;
`;

const ClearButton = styled.button`
  color: #cbd5e1;
  background: none;
  border: none;
  padding: 8px;
  border-radius: 50%;
  cursor: pointer;
  transition: all 0.2s;

  &:hover {
    color: #ef4444;
    background-color: #fef2f2;
  }
`;

const IconButton = styled.button`
  color: #94a3b8;
  background: none;
  border: none;
  padding: 8px;
  border-radius: 50%;
  cursor: pointer;
  transition: all 0.2s;

  &:hover {
    color: #0f766e;
    background-color: #f0fdfa;
  }
`;

const MessagesArea = styled.div`
  flex: 1;
  overflow-y: auto;
  padding: 20px;
  background-color: #f8fafc;
  display: flex;
  flex-direction: column;
  gap: 24px;
`;

const EmptyState = styled.div`
  display: flex;
  flex-direction: column;
  align-items: center;
  justify-content: center;
  height: 100%;
  text-align: center;
  opacity: 0.8;
  gap: 16px;
`;

const EmptyIconWrapper = styled.div`
  width: 64px;
  height: 64px;
  background-color: white;
  border-radius: 20px;
  display: flex;
  align-items: center;
  justify-content: center;
  color: #14b8a6;
  box-shadow: 0 2px 4px rgba(0,0,0,0.02);
  box-shadow: 0 2px 4px rgba(0,0,0,0.02);
  border: 1px solid #f0f9ff;
`;

const SuggestedQuestionsContainer = styled.div`
  display: flex;
  flex-direction: column;
  gap: 8px;
  width: 100%;
  margin-top: 16px;
`;

const SuggestionChip = styled.button`
  background-color: white;
  border: 1px solid #e2e8f0;
  border-radius: 12px;
  padding: 10px 14px;
  text-align: left;
  cursor: pointer;
  transition: all 0.2s;
  font-size: 13px;
  color: #475569;
  display: flex;
  align-items: center;
  gap: 8px;
  box-shadow: 0 1px 2px rgba(0,0,0,0.02);

  &:hover {
    border-color: #0f766e;
    color: #0f766e;
    background-color: #f0fdfa;
    transform: translateY(-1px);
    box-shadow: 0 4px 6px rgba(0,0,0,0.05);
  }

  svg {
    color: #94a3b8;
    transition: color 0.2s;
  }
  
  &:hover svg {
    color: #0f766e;
  }
`;

const MessageRow = styled.div`
  display: flex;
  width: 100%;
  justify-content: ${props => props.$isUser ? 'flex-end' : 'flex-start'};
`;

const MessageGroup = styled.div`
  display: flex;
  max-width: 85%;
  align-items: flex-end;
  flex-direction: ${props => props.$isUser ? 'row-reverse' : 'row'};
  gap: 8px;
`;

const BotAvatarSmall = styled.div`
  width: 24px;
  height: 24px;
  border-radius: 50%;
  background-color: #0f766e;
  display: flex;
  align-items: center;
  justify-content: center;
  color: white;
  flex-shrink: 0;
  margin-bottom: 4px;
`;

const UserAvatarSmall = styled.div`
  /* Hidden usually as per design, but defined just in case */
  width: 24px;
  height: 24px;
  border-radius: 50%;
  background-color: #4f46e5;
  display: none; 
`;

const Bubble = styled.div`
  padding: 12px 16px;
  font-size: 1em;
  line-height: 1.5;
  box-shadow: 0 1px 2px rgba(0,0,0,0.05);
  white-space: pre-wrap;
  max-width: 100%;
  
  ${props => props.$isUser ? `
    background-color: #0f766e; /* Teal-700 equivalent */
    color: white;
    border-radius: 18px 18px 2px 18px;
  ` : `
    background-color: white;
    color: #1e293b;
    border: 1px solid #e2e8f0;
    border-radius: 18px 18px 18px 2px;
  `}
`;

const ThinkingBubble = styled.div`
  background-color: white;
  border: 1px solid #e2e8f0;
  padding: 12px 16px;
  border-radius: 18px 18px 18px 2px;
  display: flex;
  align-items: center;
  gap: 6px;
  box-shadow: 0 1px 2px rgba(0,0,0,0.05);
`;

const Dot = styled.div`
  width: 4px;
  height: 4px;
  background-color: #94a3b8;
  border-radius: 50%;
  animation: bounce 1.4s infinite ease-in-out both;
  
  &:nth-child(1) { animation-delay: -0.32s; }
  &:nth-child(2) { animation-delay: -0.16s; }
  
  @keyframes bounce {
    0%, 80%, 100% { transform: scale(0); }
    40% { transform: scale(1); }
  }
`;

const InputArea = styled.div`
  padding: 16px;
  background-color: white;
  border-top: 1px solid #f1f5f9;
`;

const InputWrapper = styled.div`
  position: relative;
  /* group equivalent not needed, child focus works */
`;

const TextArea = styled.textarea`
  width: 100%;
  padding: 14px 16px;
  padding-right: 90px;
  background-color: #f8fafc;
  border: 1px solid #e2e8f0;
  border-radius: 12px;
  font-family: inherit;
  font-size: 1em;
  color: #334155;
  resize: none;
  min-height: 52px;
  outline: none;
  transition: all 0.2s;
  box-shadow: inset 0 1px 2px rgba(0,0,0,0.02);

  &:focus {
    background-color: white;
    border-color: #14b8a6;
    box-shadow: 0 0 0 3px rgba(20, 184, 166, 0.1);
  }

  &::placeholder {
    color: #94a3b8;
  }
`;

const SendButton = styled.button`
  position: absolute;
  right: 8px;
  bottom: 8px;
  padding: 8px;
  border-radius: 8px;
  border: none;
  cursor: pointer;
  transition: all 0.2s;
  display: flex;
  align-items: center;
  justify-content: center;

  ${props => props.disabled ? `
    background-color: transparent;
    color: #cbd5e1;
    cursor: not-allowed;
  ` : `
    background-color: #0f766e;
    color: white;
    box-shadow: 0 2px 4px rgba(15, 118, 110, 0.2);
    
    &:hover {
      background-color: #115e59;
      transform: translateY(-1px);
    }
    
    &:active {
      transform: translateY(0);
    }
  `}
`;

const TrashButton = styled.button`
  position: absolute;
  right: 50px;
  bottom: 8px;
  padding: 8px;
  border-radius: 8px;
  border: none;
  cursor: pointer;
  transition: all 0.2s;
  display: flex;
  align-items: center;
  justify-content: center;
  background-color: transparent;
  color: #94a3b8;

  &:hover {
    color: #ef4444;
    background-color: #fef2f2;
    transform: translateY(-1px);
  }

  &:active {
    transform: translateY(0);
  }
`;

const FooterRow = styled.div`
  margin-top: 12px;
  display: flex;
  justify-content: space-between;
  align-items: center;
  padding: 0 4px;
`;

const SelectWrapper = styled.div`
  position: relative;
  display: flex;
  align-items: center;
  
  &:hover svg {
    color: #0f766e;
  }
`;

const ModelSelect = styled.select`
  appearance: none;
  background-color: transparent;
  border: none;
  padding: 4px 24px 4px 0;
  font-size: 12px;
  font-weight: 600;
  color: #64748b;
  cursor: pointer;
  outline: none;
  transition: color 0.2s;

  &:hover {
    color: #0f766e;
  }
`;

const IconWrapper = styled.div`
  position: absolute;
  right: 0;
  pointer-events: none;
  color: #94a3b8;
  transition: color 0.2s;
`;

const Disclaimer = styled.span`
  font-size: 10px;
  color: #cbd5e1;
  font-weight: 500;
  font-size: 10px;
  color: #cbd5e1;
  font-weight: 500;
`;

const ClearChatLink = styled.button`
  background: none;
  border: none;
  color: #94a3b8;
  font-size: 11px;
  cursor: pointer;
  margin-top: 8px;
  align-self: flex-end;
  text-decoration: underline;
  padding: 4px 8px;
  transition: color 0.2s;
  
  &:hover {
    color: #ef4444;
  }
`;

// --- Component ---

const Chatbot = () => {
  const [messages, setMessages] = useState([]);
  const [input, setInput] = useState('');
  const [isLoading, setIsLoading] = useState(false);
  const [selectedModel, setSelectedModel] = useState('gpt'); // 'gpt' or 'gemini'
  const [isMinimized, setIsMinimized] = useState(() => {
    const savedState = localStorage.getItem('chatbot_minimized');
    return savedState === 'true';
  });
  const [isDocked, setIsDocked] = useState(false);
  const [copiedIndex, setCopiedIndex] = useState(null);
  const messagesEndRef = useRef(null);
  
  const PRESET_QUESTIONS = [
    { title: "What cell type does each cluster most likely represent using layer: MAGIC?", "prompt": "Identify cell types of Aorta cells from mouse using the following markers separately for each row. Some can be a mixture of multiple cell types.\n GeneList’:\nCluster 0: Rarres2, Mmp23, Col3a1, Col6a1, C1s, Gas1, Aebp1, Olfml2b, Col1a2, Htra1, Pcolce, Col5a2, Lhfp, Ddr2, Serping1, Chpf, C1ra, Col6a2, Dcn, Timp2, Ptgis, Nupr1, Ccdc80, Sod3, Col1a1, Igfbp6, Plxdc2, Lgals1, C3, Fndc1, Loxl1, Cpxm1, Cfh, C1qtnf2, C4b, Rcn3, Cercam, Mgp, Serpinf1, Tmem119, Bmp1, Srpx, Cp, Fn1, Stbd1, Olfml3, Loxl3, Chrdl1, Scd1, Eln\nCluster 1: Aqp7, Rbp7, 8430408G22Rik, Magix, C1qtnf9, Slc26a10, Car4, Pparg, Timp4, Meox2, Cd36, Fam70b, Pbld1, Itga1, Kif26a, Nepn, Cxcl12, Prdm16, Cyp1a1, Hspa12b, St6galnac2, P2ry2, Ptprr, Pkp4, Cav1, Fabp5, Cdh13, Gfod1, Apold1, AW112010, Jam2, Hey1, Cxcl9, Fam176a, Magi1, Cldn5, Gpr160, Palmd, Eepd1, Myzap, Sdpr, Aqp1, Sox17, Fmo1, Rasd1, Lipe, Rnf125, Snx32, N4bp3, Tspan7\nCluster 2: Esm1, Ehd3, Lrg1, Gpr97, Dkk2, Plvap, Igfbp3, Mcam, Ntn4, Lrrc3b, Cd1d1, Pydc3, Ica1, Lama3, Dchs1, Hapln1, Plaur, Fam40a, Klhl2, Rcsd1, Ccdc67, Fam171a1, Hps6, Oasl1, Tnfaip8l1, I830012O16Rik, Fam174b, Dcbld1, Ankrd29, 4933407C03Rik, Yes1, Cldn15, Dffb, Ppm1f, Efna1, Kctd12b, Mcat, Ada, Acta1, 2010321M09Rik, Scarf1, Tbc1d20, Mafb, Spry4, Card6, Mcm3, Tbc1d9b, Stab1, Arrdc2, Peg3\nCluster 3: Cst6, Tbx21, 4933427D14Rik, Dzip1, Fggy, Rnmtl1, Gpc1, Epb4.1l4a, 5730590G19Rik, Aass, Cenpf, Fhl2, 6330403A02Rik, S100pbp, Cd8b1, Acta2, Ccl22, Gzma, Ucp1, C3ar1, Cd96, Gm525, Hspb6, Adck3, Ccrl1, Gpd1, Ubash3a, Cd3g, Fbxw10, Palb2, Myl7, Muc5b, Akr1c12, Myoc, Tcap, Car3, Ms4a4b, Myom1, Slc43a2, Ccna2, Cenpe, Ticam1, Dscc1, 9930013L23Rik, Cd247, Ccl5, Fbxl2, Tnnt2, Has1, Myh11\nCluster 4: Tyrobp, Fcgr2b, Cd37, Cd83, Ctss, Ptpn6, Ctsc, Il10ra, Cfp, Prkcb, Cd68, Adrb2, Cyth4, H2-Aa, Epsti1, Aif1, H2-Eb1, Arhgef6, Lpxn, H2-DMa, Ly86, Laptm5, Pld4, Cd74, Sema4a, Ptprc, H2-Ab1, Ms4a6c, Bcl2a1b, Rgs14, Cd300a, Rassf4, Gm11428, Lyz2, Sh3kbp1, Lat2, Rac2, P2ry6, Itgb2, Zc3h12d, Rasal3, Hck, Ncf2, Lcp2, Plbd1, Ccl6, Lcp1, Coro1a, Tep1, Sgpl1\nCluster 5: Cyp2f2, Wfdc2, Tspan1, Krt5, Reg3g, Krt15, Bpifa1, Aqp3, Ckmt1, Ifitm1, 5330417C22Rik, Krt8, Anxa8, Ces1d, Aldh3a1, Elf3, Scgb1a1, Krt18, Fmo3, Plxnb1, Krt7, Emb, Car5b, Bpifb1, Niacr1, Aqp5, Atp1b1, Aox3, Gprc5a, Adh7, Sdc1, Ehf, Irx2, Cbr2, Ocln, Plcd3, Zbtb49, Krt17, Slc22a23, Cldn8, Acpp, Krt19, Pdgfc, Prodh, F3, Epcam, Fam187b, Tjp3, Ccdc51, Ccdc3\nCluster 6: Zfp39, Rnmtl1, Samd10, Nat9, Zfp64, 4933439C10Rik, Gbp3, Zfp747, Fabp4, 9930014A18Rik, Rnf152, Cldn5, Tbx21, Acacb, Cav1, 5730590G19Rik, Tcf7, Hist1h2be, Tm6sf2, Tox, Prodh, Cd36, Ucp1, Nqo1, Pecam1, Epb4.1l4a, Gzma, Aqp1, Cytl1, Sdcbp2, Pkn3, Slc52a2, Sdpr, Enpp6, Lyve1, Chaf1b, Akap2, Mtm1, Isl1, Fam70b, Car3, Selp, Apol9a, Clu, Myf5, Pde9a, Gata3, Mgp, Hey2, Fmo1\nCluster 7: Hba-a1, Beta-s, Apol11b, Hbb-b1, Alas2, Slc25a37, Fam46c, Isg20, Myom1, Snca, Fech, Epb4.1, Bpgm, Mkrn1, Ube2l6, 2810453I06Rik, Bco2, Pira2, Srl, Nusap1, Foxred2, Kcp, Clec4a4, Tspo2, Cd4, Melk, Col2a1, Cdc25b, Gypa, 5730469M10Rik, Eomes, Tlr12, Trim10, Car2, Ube2c, Eaf2, Mdga1, Adipoq, Cdr2, Sec1, Aldh1b1, Ms4a4b, Atad5, Scin, Coch, Col9a2, Cenpf, Traf3ip3, St8sia1, Ccndbp1\nProvide the most likely cell type for each cluster based on these marker genes and show your reasoning." },
  ];
  console.log("PRESET_QUESTIONS: ", PRESET_QUESTIONS);
  useEffect(() => {
    localStorage.setItem('chatbot_minimized', isMinimized);
  }, [isMinimized]);

  // Resize state
  const [size, setSize] = useState({ width: 380, height: 600 });
  const isResizing = useRef(false);
  const startPos = useRef({ x: 0, y: 0 });
  const startSize = useRef({ w: 0, h: 0 });

  const scrollToBottom = () => {
    messagesEndRef.current?.scrollIntoView({ behavior: "smooth" });
  };

  useEffect(() => {
    scrollToBottom();
  }, [messages]);

  useEffect(() => {
    const handleMouseMove = (e) => {
      if (!isResizing.current) return;
      const deltaX = e.clientX - startPos.current.x;
      const deltaY = e.clientY - startPos.current.y;

      if (Number.isFinite(startSize.current.w - deltaX) && Number.isFinite(startSize.current.h - deltaY)) {
        setSize({
          width: Math.max(300, startSize.current.w - deltaX),
          height: Math.max(400, startSize.current.h - deltaY)
        });
      }
    };

    const handleMouseUp = () => {
      isResizing.current = false;
      document.body.style.cursor = 'default';
      document.body.style.userSelect = 'auto';
    };

    document.addEventListener('mousemove', handleMouseMove);
    document.addEventListener('mouseup', handleMouseUp);

    return () => {
      document.removeEventListener('mousemove', handleMouseMove);
      document.removeEventListener('mouseup', handleMouseUp);
    };
  }, []);

  const startResize = (e) => {
    e.preventDefault();
    isResizing.current = true;
    startPos.current = { x: e.clientX, y: e.clientY };
    startSize.current = { w: size.width, h: size.height };
    document.body.style.cursor = 'nw-resize';
    document.body.style.userSelect = 'none';
  };

  const handleSend = async () => {
    if (!input.trim()) return;

    const userMessage = { role: 'user', content: input };
    setMessages(prev => [...prev, userMessage]);
    setInput('');
    setIsLoading(true);

    const targetModel = selectedModel; // Capture current model

    try {
      const response = await axios.post(`${NODE_API_URL}/api/chat`, {
        message: userMessage.content,
        model: targetModel
      });

      const botMessage = { role: 'assistant', content: response.data.reply };
      setMessages(prev => [...prev, botMessage]);
    } catch (error) {
      console.error("Chat error:", error);
      let messageContent = "Sorry, something went wrong. Please check your API keys in oscb-node/.env.";

      if (error.response && error.response.status === 429) {
        messageContent = "You have exceeded the API quota (Rate Limit). Please wait a moment before trying again.";
      }

      const errorMessage = { role: 'assistant', content: messageContent };
      setMessages(prev => [...prev, errorMessage]);
    } finally {
      setIsLoading(false);
    }
  };

  const handleClear = () => {
    setMessages([]);
  };

  const handleKeyDown = (e) => {
    if (e.key === 'Enter' && !e.shiftKey) {
      e.preventDefault();
      handleSend();
    }
  };

  if (isMinimized) {
    return (
      <MinimizedButton onClick={() => setIsMinimized(false)} title="Open AI Assistant">
        <img src={SingleCellLogo} alt="AI" style={{ width: '100%', height: '100%', borderRadius: '50%' }} />
      </MinimizedButton>
    );
  }

  const handleCopy = (index) => {
    setCopiedIndex(index);
    setTimeout(() => setCopiedIndex(null), 2000); // Reset after 2 seconds
  };

  let codeBlockIndex = -1;

  return (
    <Container
      style={{
        width: isDocked ? '600px' : `${size.width}px`,
        height: isDocked ? '500px' : `${size.height}px`,
        right: isDocked ? '24px' : '24px',
        bottom: isDocked ? '0px' : '24px',
        borderBottomRightRadius: isDocked ? '0' : '16px',
        borderBottomLeftRadius: isDocked ? '0' : '16px',
        fontSize: size.width > 600 ? '18px' : size.width > 450 ? '16px' : '14px' // Enhanced responsive font size
      }}
    >
      {!isDocked && <ResizeHandle onMouseDown={startResize} title="Drag to resize" />}
      <Header>
        <HeaderTitleGroup>
          <AvatarCircle style={{ background: 'transparent', boxShadow: 'none' }}>
            <img src={SingleCellLogo} alt="AI" style={{ width: '100%', height: '100%', borderRadius: '50%' }} />
          </AvatarCircle>
          <div>
            <TitleText>AI Assistant</TitleText>
            <StatusIndicator>
              <StatusDot />
              <StatusText>Online</StatusText>
            </StatusIndicator>
          </div>
        </HeaderTitleGroup>

        <HeaderActions>
          <IconButton onClick={() => setIsMinimized(true)} title="Minimize">
            <FontAwesomeIcon icon={faMinus} size="sm" />
          </IconButton>
        </HeaderActions>
      </Header>

      <MessagesArea>
        {messages.length === 0 && PRESET_QUESTIONS && (
          <EmptyState>
            <EmptyIconWrapper>
              <FontAwesomeIcon icon={faMagic} size="lg" />
            </EmptyIconWrapper>
            <div style={{ display: 'flex', flexDirection: 'column', gap: '4px' }}>
              <p style={{ margin: 0, fontSize: '16px', fontWeight: 600, color: '#334155' }}>How can I help you?</p>
              <p style={{ margin: 0, fontSize: '12px', color: '#94a3b8', maxWidth: '220px', lineHeight: '1.5' }}>
                I can answer questions about single-cell sequencing analysis. This AI assistant may occasionally generate incorrect or misleading information. We are not responsible for any decisions made based on the generated content. Please verify critical information independently.
              </p>
            </div>

            <SuggestedQuestionsContainer>
              {Array.isArray(PRESET_QUESTIONS) && PRESET_QUESTIONS.map((q, idx) => (
                <SuggestionChip key={idx} onClick={() => setInput(q.prompt)}>
                  <FontAwesomeIcon icon={faMagic} size="xs" />
                  {q.title}
                </SuggestionChip>
                  )
                )
              }
            </SuggestedQuestionsContainer>
          </EmptyState>
        )}

        {messages.map((msg, index) => (
          <MessageRow key={index} $isUser={msg.role === 'user'}>
            <MessageGroup $isUser={msg.role === 'user'}>
              {/* Bot Avatar */}
              {msg.role !== 'user' && (
                <BotAvatarSmall>
                  <FontAwesomeIcon icon={faRobot} size="xs" />
                </BotAvatarSmall>
              )}

              {/* Bubble */}
              <Bubble $isUser={msg.role === 'user'}>
                <ReactMarkdown
                  remarkPlugins={[remarkGfm]}
                  rehypePlugins={[rehypeRaw, rehypeGithubAlerts]}
                  children={msg.content}
                  components={{
                    code(props) {
                      const { children, inline, className, node, ...rest } = props;
                      const match = /language-(\w+)/.exec(className || '');
                      const codeText = String(children).replace(/\n$/, '');
                      if (!inline && match) {
                        codeBlockIndex++;

                        const currentIndex = codeBlockIndex;
                        const codeText = String(children).replace(/\n$/, '');

                        return (
                          <div style={{ position: 'relative' }}>
                            <SyntaxHighlighter
                              {...rest}
                              PreTag="div"
                              children={codeText}
                              language={match[1]}
                              style={dark}
                            />
                            <CopyToClipboard text={codeText} onCopy={() => handleCopy(currentIndex)}>
                              <button style={{
                                position: 'absolute',
                                top: '5px',
                                right: '5px',
                                background: '#333',
                                color: '#fff',
                                border: 'none',
                                borderRadius: '4px',
                                cursor: 'pointer',
                                padding: '5px 10px',
                              }}>{copiedIndex === currentIndex ? 'Copied!' : 'Copy'}</button>
                            </CopyToClipboard>
                          </div>
                        );
                      }

                      // Fallback for inline code or unknown language
                      return (
                        <code {...rest} className={className}>
                          {children}
                        </code>
                      );
                    }
                  }}
                />
              </Bubble>

              {/* User Avatar (Hidden) */}
              {msg.role === 'user' && (
                <UserAvatarSmall>
                  <FontAwesomeIcon icon={faUser} size="xs" />
                </UserAvatarSmall>
              )}
            </MessageGroup>
          </MessageRow>
        ))}

        {isLoading && (
          <MessageRow $isUser={false}>
            <MessageGroup $isUser={false}>
              <BotAvatarSmall>
                <FontAwesomeIcon icon={faRobot} size="xs" />
              </BotAvatarSmall>
              <ThinkingBubble>
                <span style={{ fontSize: '12px', color: '#94a3b8', fontWeight: 500 }}>Thinking</span>
                <div style={{ display: 'flex', gap: '2px', marginLeft: '2px' }}>
                  <Dot />
                  <Dot />
                  <Dot />
                </div>
              </ThinkingBubble>
            </MessageGroup>
          </MessageRow>
        )}
        <div ref={messagesEndRef} />
      </MessagesArea>

      <InputArea>
        <InputWrapper>
          <TextArea
            value={input}
            onChange={(e) => setInput(e.target.value)}
            onKeyDown={handleKeyDown}
            placeholder="Type your message..."
            rows="1"
          />
          <TrashButton
            onClick={handleClear}
            title="Clear chat history"
            style={{ display: messages.length > 0 ? 'flex' : 'none' }}
          >
            <FontAwesomeIcon icon={faTrash} size="sm" />
          </TrashButton>
          <SendButton
            onClick={handleSend}
            disabled={isLoading || !input.trim()}
          >
            <FontAwesomeIcon icon={faPaperPlane} size="sm" />
          </SendButton>
        </InputWrapper>

        <FooterRow>
          <SelectWrapper>
            <ModelSelect
              value={selectedModel}
              onChange={(e) => setSelectedModel(e.target.value)}
            >
              <option value="gpt">GPT-4o (OpenAI)</option>
              <option value="gemini">Gemini Flash (Latest) (Google)</option>
            </ModelSelect>
            <IconWrapper>
              <FontAwesomeIcon icon={faRotateRight} rotation={90} size="xs" />
            </IconWrapper>
          </SelectWrapper>
          <Disclaimer>Powered by AI</Disclaimer>
        </FooterRow>
      </InputArea>
    </Container >
  );
};

export default Chatbot;
