import { PageContainer } from '@ant-design/pro-components';
import React, { useState, useEffect, useRef } from 'react';
import { Card, Input, Button, Typography } from 'antd';
import './predict.css'; // 导入 CSS 文件

const Predict: React.FC = () => {
  const [smiles, setSmiles] = useState('');
  const [output, setOutput] = useState('预测结果将在此处显示');
  const viewerRef = useRef<HTMLDivElement | null>(null);
  
  // 顶层 useEffect，用于动态加载 3Dmol.js
  useEffect(() => {
    const script = document.createElement('script');
    script.src = `../../../public/scripts/3Dmol.js`;
    script.onload = () => {
      console.log('3Dmol.js loaded');
    };
    document.body.appendChild(script);

    // 清理脚本
    return () => {
      document.body.removeChild(script);
    };
  }, []); // 空依赖数组，确保仅加载一次

  // 定义 loadMolecule 函数
  const loadMolecule = () => {
    // 检查 3Dmol.js 是否加载完成
    if (!(window as any).$3Dmol) {
      console.log("3Dmol.js 未加载，无法创建查看器实例");
      return;
    }
  
    // 检查 viewerRef 是否已经初始化
    if (!viewerRef.current) {
      console.log("viewerRef 未初始化");
      return;
    }
  
    console.log("3Dmol.js 已加载，viewerRef 已初始化");
  
    // 创建 3Dmol 查看器实例
    const viewer = (window as any).$3Dmol.createViewer(viewerRef.current, {
      backgroundColor: 'white',
    });
    console.log("3Dmol viewer 实例已创建");
  
    // 发送请求以获取分子数据
    fetch(`http://127.0.0.1:5273/generate_mol?smiles=${encodeURIComponent(smiles)}`)
      .then(response => {
        if (!response.ok) {
          throw new Error(`网络请求失败，状态码: ${response.status}`);
        }
        console.log("网络请求成功，正在读取分子数据...");
        return response.text();
      })
      .then(molData => {
        console.log("分子数据加载成功:", molData);
  
        // 添加分子模型到 viewer
        viewer.addModel(molData, 'mol');
        viewer.setStyle({}, { stick: {} });
        viewer.zoomTo();
        viewer.render();
        console.log("分子模型已渲染");
      })
      .catch(error => {
        console.error("加载分子模型时出错:", error);
      });
  };
    const handleButtonClick = () => {
      loadMolecule();
      handlePredict();
    };
  const handlePredict = () => {
    setOutput('正在预测属性...');
    
    fetch('http://localhost:9000/predict', {
      method: 'POST',
      mode: 'cors',
      headers: {
        'Content-Type': 'application/json'
      },
      body: JSON.stringify({ smiles })
    })
    
      .then(response => {
        if (!response.ok) {
          throw new Error("网络响应错误，状态码: " + response.status);
        }
        return response.text();  // 将响应作为纯文本读取
      })
      .then(text => {
        try {
          const data = JSON.parse(text);  // 尝试解析 JSON 数据
          if (data.error) {
            setOutput("预测失败: " + data.error);
          } else {
            setOutput(`
              预测结果:
              分子量: ${data.data.分子量 || 'N/A'} g/mol
              体积: ${data.data.体积 || 'N/A'} Å³
              密度: ${data.data.密度 || 'N/A'} g/cm³
              氢键受体数: ${data.data.氢键受体数 || 'N/A'}
              氢键供体数: ${data.data.氢键供体数 || 'N/A'}
              可旋转键数: ${data.data.可旋转键数 || 'N/A'}
              环数: ${data.data.环数 || 'N/A'}
              杂原子数：${data.data.杂原子数 || 'N/A'}
              TPSA: ${data.data.TPSA || 'N/A'}
              logP: ${data.data.logP || 'N/A'}
              logS: ${data.data.logS || 'N/A'}
              毒性: ${data.data.toxicity || 'N/A'}
            `);
          }
        } catch (error) {
          console.error("Failed to parse JSON:", error);
          setOutput("预测失败，服务器返回了非 JSON 响应：" + text);
        }
      })
      .catch(error => {
        console.error("Error during fetch:", error);
        setOutput("预测失败，请重试。详细错误信息：" + error.message);
      });
  };
  

  return (
    <PageContainer>
      <Card>
        <h2>分子属性预测页面</h2>
        <div className="input-section">
          <label htmlFor="smiles">SMILES:</label>
          <Input
            id="smiles"
            value={smiles}
            onChange={(e) => setSmiles(e.target.value)}
            placeholder="请输入 SMILES 字符串"
          />
          <Button type="primary" onClick={handleButtonClick} >确认目标化合物</Button>
        </div>
         {/* 空白容器，展示 3D 分子模型 */}
         <div style={{ margin: '20px 0' }}>
          <Card title="3D 分子模型" bordered={false}>
            <div
              id="viewer"
              ref={viewerRef}
              style={{ width: '100%', height: '400px', position: 'relative' }}
            ></div>
          </Card>
        </div>
        <h3>预测结果</h3>
        <div className="output-box">
          <pre>{output}</pre>
        </div>
      </Card>
    </PageContainer>
  );
};

export default Predict;
