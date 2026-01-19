import React, { useState, useRef, useEffect } from 'react';
import axios from 'axios';
import { FontAwesomeIcon } from '@fortawesome/react-fontawesome';
import { faPaperPlane, faTrash, faUser, faRobot, faRotateRight, faMagic, faMinus, faExpand } from '@fortawesome/free-solid-svg-icons';
import { NODE_API_URL } from '../../constants/declarations';
import styled from 'styled-components';

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
  background: linear-gradient(135deg, #0f766e 0%, #0d9488 100%);
  color: white;
  border: none;
  box-shadow: 0 4px 12px rgba(13, 148, 136, 0.4);
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
  border: 1px solid #f0f9ff;
`;

const MessageRow = styled.div`
  display: flex;
  width: 100%;
  justify-content: ${props => props.$isUser ? 'flex-end' : 'flex-start'};
`;

const MessageGroup = styled.div`
  display: flex;
  /* max-width: 85%; */
  /* Using inline styles for responsive width if needed, but 85% is good default */
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
  font-size: 14px;
  line-height: 1.5;
  box-shadow: 0 1px 2px rgba(0,0,0,0.05);
  white-space: pre-wrap;
  max-width: 280px; /* Constrain width */
  
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
  padding-right: 48px;
  background-color: #f8fafc;
  border: 1px solid #e2e8f0;
  border-radius: 12px;
  font-family: inherit;
  font-size: 14px;
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
`;

// --- Component ---

const Chatbot = () => {
  const [messages, setMessages] = useState([]);
  const [input, setInput] = useState('');
  const [isLoading, setIsLoading] = useState(false);
  const [selectedModel, setSelectedModel] = useState('gpt'); // 'gpt' or 'gemini'
  const [isMinimized, setIsMinimized] = useState(false);
  const messagesEndRef = useRef(null);

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

      setSize({
        width: Math.max(300, startSize.current.w - deltaX),
        height: Math.max(400, startSize.current.h - deltaY)
      });
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

    try {
      const response = await axios.post(`${NODE_API_URL}/api/chat`, {
        message: userMessage.content,
        model: selectedModel
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
        <FontAwesomeIcon icon={faRobot} size="lg" />
      </MinimizedButton>
    );
  }

  return (
    <Container style={{ width: `${size.width}px`, height: `${size.height}px` }}>
      <ResizeHandle onMouseDown={startResize} title="Drag to resize" />
      <Header>
        <HeaderTitleGroup>
          <AvatarCircle>
            <FontAwesomeIcon icon={faRobot} size="sm" />
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
          <ClearButton onClick={handleClear} title="Clear Conversation">
            <FontAwesomeIcon icon={faTrash} size="sm" />
          </ClearButton>
        </HeaderActions>
      </Header>

      <MessagesArea>
        {messages.length === 0 && (
          <EmptyState>
            <EmptyIconWrapper>
              <FontAwesomeIcon icon={faMagic} size="lg" />
            </EmptyIconWrapper>
            <div style={{ display: 'flex', flexDirection: 'column', gap: '4px' }}>
              <p style={{ margin: 0, fontSize: '16px', fontWeight: 600, color: '#334155' }}>How can I help you?</p>
              <p style={{ margin: 0, fontSize: '12px', color: '#94a3b8', maxWidth: '220px', lineHeight: '1.5' }}>
                I can answer questions about single-cell analysis, datasets, and RNA sequencing.
              </p>
            </div>
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
                {msg.content}
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
              <option value="gemini">Gemini Pro (Google)</option>
            </ModelSelect>
            <IconWrapper>
              <FontAwesomeIcon icon={faRotateRight} rotation={90} size="xs" />
            </IconWrapper>
          </SelectWrapper>
          <Disclaimer>Powered by AI</Disclaimer>
        </FooterRow>
      </InputArea>
    </Container>
  );
};

export default Chatbot;
