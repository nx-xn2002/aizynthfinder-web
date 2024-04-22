import Chat, { Bubble, useMessages } from '@chatui/core';
import '@chatui/core/dist/index.css';
import axios from 'axios';
import {useState} from 'react';
import '@chatui/core/es/styles/index.less';
import { InputNumber } from 'antd';
import { Button, Flex } from 'antd';
// 引入样式
import '@chatui/core/dist/index.css';
import { Input } from 'antd';
const initialMessages = [
  {
    type: 'text',
    content: { text: '您好，我是智能助理，您的贴心小助手~' },
    user: { avatar: '/logo.svg' },
  },
  // {
  //   type: 'image',
  //   content: {
  //     picUrl: '//img.alicdn.com/tfs/TB1p_nirYr1gK0jSZR0XXbP8XXa-300-300.png',
  //   },
  // },
];

// 默认快捷短语，可选
const defaultQuickReplies = [
  // {
  //   icon: 'message',
  //   name: '联系人工服务',
  //   isNew: true,
  //   isHighlight: true,
  // },
  //   {
  //     name: '短语1',
  //     isNew: true,
  //   },
  //   {
  //     name: '短语2',
  //     isHighlight: true,
  //   },
  //   {
  //     name: '短语3',
  //   },
];

const ChatBot = () => {
  // 消息列表
  const { messages, appendMsg, setTyping } = useMessages(initialMessages);
  const [inputValue, setInputValue] = useState("CC(=O)OC1=CC=C(C=C1)C(=O)OCC");
  const [temperature, setTemperature] = useState("0.1");
  const [top_p, setTopp] = useState("0.75");
  const [top_k, setTopk] = useState("40");
  const [num_beams, setNumBeams] = useState("4");
  const [repetition_penalty, setRepetitionPenalty] = useState("1");
  const [max_new_tokens, setMaxNewTokens] = useState("128");
  const handleTemperatureChange = (value: string)  => {
    setTemperature(value);
  };
  const handleToppChange = (value: string)   => {
    console.log('top_p_change:', value);
    setTopp(value);
  };
  const handleTopkChange = (value: string)  => {
    setTopk(value);
  };
  const handleNumBeamsChange = (value: string)  =>  {
    setNumBeams(value);
  };
  const handleRepetitionPenaltyChange = (value: string)  => {
    setRepetitionPenalty(value);
  };
  const handleMaxNewTokensChange = (value: string)  => {
    setMaxNewTokens(value);
  };
  const handleReset = () => {
    setTemperature("0.1"); // 重置为初始值 "1"
    setTopp("0.75");
    setTopk("40");
    setNumBeams("4");
    setRepetitionPenalty("1");
    setMaxNewTokens("128");
  };
  // 发送回调
  function handleSend(type, val) {
    console.log('type-val:', type, '-', val);
    // console.log("Sending message:", inputValue);
    if (type === 'text' && val.trim()) {
      // TODO: 发送请求
      appendMsg({
        type: 'text',
        content: { text: val },
        position: 'right',
      });
      if(val === '联系人工服务') {
        appendMsg({
          type: 'text',
          content: {text: '请拨打电话：123456789'},
          user: {avatar: '/logo.svg'},
        });
      } else {
        setTyping(true);
        axios
          .post(
            'http://localhost:9000/associate/evaluate',
            {
              input_text: inputValue,
              temperature: temperature,
              top_p: top_p,
              top_k: top_k,
              num_beams: num_beams,
              repetition_penalty: repetition_penalty,
              max_new_tokens: max_new_tokens,
              instruction: val,
            },
            {
              headers: {
                'Content-Type': 'application/json',
              },
            },
          )
          .then((res) => {
            const response = res.data.data.result;
            appendMsg({
              type: 'text',
              content: { text: response },
              user: { avatar: '/logo.svg' },
            });
          })
          .catch((err) => console.log(err));
      }
    }
  }

  // 快捷短语回调，可根据 item 数据做出不同的操作，这里以发送文本消息为例
  function handleQuickReplyClick(item) {
    handleSend('text', item.name);
  }

  function renderMessageContent(msg) {
    const { type, content } = msg;
    // 根据消息类型来渲染
    switch (type) {
      case 'text':
        return <Bubble content={content.text} />;
      case 'image':
        return (
          <Bubble type="image">
            <img src={content.picUrl} alt="" />
          </Bubble>
        );
      default:
        return null;
    }
  }

   const handleChange = (e) => {
    setInputValue(e.target.value);
  };

  return (
    <div style={{ display: 'flex', flexDirection: 'column', marginBottom: '20px' }}>
      <Chat
          navbar={{ title: '药研助手' }}
          messages={messages}
          renderMessageContent={renderMessageContent}
          quickReplies={defaultQuickReplies}
          onQuickReplyClick={handleQuickReplyClick}
          onSend={handleSend}
      />
      <div className="cool-input" style={{ marginBottom: '10px', display: 'flex', alignItems: 'center' }}>
          <label htmlFor="molecularFormula" style={{ color: 'black', marginRight: '10px',fontWeight: 'bold' }}>分子式</label>
          <Input
              id="molecularFormula"
              placeholder="Enter molecular formula"
              value={inputValue}
              onChange={handleChange}
              style={{ width: '300px' }}
          />
      </div>
      <div className="cool-input" style={{ marginBottom: '10px',display: 'flex', alignItems: 'center' }}>
          <label  style={{ color: 'black', marginRight: '10px',fontWeight: 'bold' }}>Temperature</label>
          <InputNumber<string>
              id="temperature"
              style={{ width: '200px' }}
              value={temperature}
              min="0"
              max="10"
              step="0.1"
              onChange={handleTemperatureChange}
              stringMode
          />
      </div>
      <div className="cool-input" style={{ marginBottom: '10px',display: 'flex', alignItems: 'center' }}>
          <label  style={{ color: 'black', marginRight: '10px',fontWeight: 'bold' }}>top_p</label>
          <InputNumber<string>
              id="top_p"
              style={{ width: '200px' }}
              defaultValue={top_p}
              value={top_p}
              min="0"
              max="100"
              step="0.01"
              onChange={handleToppChange}
              stringMode
          />
      </div>
      <div className="cool-input" style={{ marginBottom: '10px',display: 'flex', alignItems: 'center' }}>
          <label htmlFor="temperature" style={{ color: 'black', marginRight: '10px',fontWeight: 'bold' }}>top_k</label>
          <InputNumber<string>
              style={{ width: '200px' }}
              defaultValue={top_k}
              value={top_k}
              min="0"
              max="100"
              step="0.01"
              onChange={handleTopkChange}
              stringMode
          />
      </div>
       <div className="cool-input" style={{ marginBottom: '10px',display: 'flex', alignItems: 'center' }}>
          <label  style={{ color: 'black', marginRight: '10px',fontWeight: 'bold' }}>num_beams</label>
          <InputNumber<string>
              style={{ width: '200px' }}
              defaultValue={num_beams}
              value={num_beams}
              min="0"
              max="10"
              step="1"
              onChange={handleNumBeamsChange}
              stringMode
          />
      </div>
      <div className="cool-input" style={{ marginBottom: '10px',display: 'flex', alignItems: 'center' }}>
          <label  style={{ color: 'black', marginRight: '10px',fontWeight: 'bold' }}>repetition_penalty</label>
          <InputNumber<string>
              style={{ width: '200px' }}
              defaultValue={repetition_penalty}
              value={repetition_penalty}
              min="0"
              max="10"
              step="1"
              onChange={handleRepetitionPenaltyChange}
              stringMode
          />
      </div>
       <div className="cool-input" style={{ marginBottom: '10px',display: 'flex', alignItems: 'center' }}>
          <label  style={{ color: 'black', marginRight: '10px',fontWeight: 'bold' }}>max_new_tokens</label>
          <InputNumber<string>
              style={{ width: '200px' }}
              defaultValue={max_new_tokens}
              value={max_new_tokens}
              min="0"
              max="2048"
              step="1"
              onChange={handleMaxNewTokensChange}
              stringMode
          />
      </div>
      <Flex gap="small" wrap="wrap">
        <Button type="primary" onClick={handleReset}>重置</Button>
      </Flex>
    </div>
  );
};

export default ChatBot;
