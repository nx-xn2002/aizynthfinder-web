import {PageContainer,} from '@ant-design/pro-components';
import React, {useEffect, useState} from 'react';
import {Button, Input, Table, message, Typography, Select} from 'antd';
import {IngredientQueryRequest, listTcmByPageUsingPost} from "@/services/tcm/TcmController";
import {ColumnsType} from "antd/es/table";
import { history } from 'umi';
const { Text } = Typography;
const columns: ColumnsType<API.TcmResponse> = [
  {
    title: '中文名称',
    dataIndex: 'chineseName',
    key: 'chineseName',
    render: (text) => <Text ellipsis>{text}</Text>
  },
  {
    title: '英文名称',
    dataIndex: 'englishName',
    key: 'englishName',
  },
  {
    title: '拼音名称',
    dataIndex: 'pinyinName',
    key: 'pinyinName',
  },
  {
    title: '使用部位',
    dataIndex: 'usePart',
    key: 'usePart',
  },
  {
    title: '功效',
    dataIndex: 'effect',
    key: 'effect',
  },
  {
    title: '性味',
    dataIndex: 'propertyFlavor',
    key: 'propertyFlavor',
  },
  {
    title: '作用区域',
    dataIndex: 'channelTropism',
    key: 'channelTropism',
  },
  {
    title: '功能主治',
    dataIndex: 'indication',
    key: 'indication',
  },
  {
    title: '操作',
    key: 'action',
    render: (txt, record, index) =>{
      const handleClick = async () => {
        //console.log('按钮点击，当前行数据对象：', record);
        const res = await IngredientQueryRequest(record);
        //console.log('响应', res.data.list);
        history.push('/generate', {ingredients : res.data.list,id:record.id});
      };
      return (
        <div>
          <Button type="primary" size="small" onClick={handleClick}>生成</Button>&nbsp;
        </div>
      )
    }
  },
];
const Tcm: React.FC = () => {
  const [loading, setLoading] = useState(false);
  const { Option } = Select;
  const [searchType, setSearchType] = useState('chinese_name');
  const [searchValue, setSearchValue] = useState('');

  const initSearchParams = {
    current: 1,
    pageSize: 10,
    chinese_name: '',
    english_name: '',
    pinyin_name: '',
  };
  // 中药查询参数
  const [searchParams] = useState<API.TcmQueryRequest>({...initSearchParams});
  const [tcmList, setTcmListList] = useState<API.TcmResponse[]>();
  const [total, setTotal] = useState<number>(0);

  //搜索框
  const handleSearch = async () => {
    setLoading(true);
    try {
      // 构建请求参数
      const requestData = {
        [searchType]: searchValue
      };
      // 发送post请求
      const res = await listTcmByPageUsingPost(requestData); // 假设listTcmByPageUsingPost是发送post请求的函数
      // 请求成功并得到响应结果
      setLoading(false);
      setTcmListList(res.data.list ?? []);
      setTotal(res.data.total ?? 0);
      message.log('请求成功');
    } catch (error) {
      // 请求失败
      setLoading(false);
      console.error('请求失败:', error);
    }
  };

  const loadData = async (page: any, pageSize: any) => {
    try {
      searchParams.pageSize = pageSize ?? 10;
      searchParams.current = page ?? 1;
      // 检查是否有搜索值
      if (searchType==='chinese_name') {
        // 如果有搜索值，则添加到搜索参数中
        searchParams.chinese_name = searchValue;
      } else if(searchType==='english_name') {
        // 如果没有搜索值，则移除搜索参数中的 searchValue
        searchParams.english_name = searchValue;
      }else if(searchType==='pinyin_name') {
        // 如果没有搜索值，则移除搜索参数中的 searchValue
        searchParams.pinyin_name = searchValue;
      }
      const res = await listTcmByPageUsingPost(searchParams);
      if (res.data) {
        setTcmListList(res.data.list ?? []);
        setTotal(res.data.total ?? 0);
      } else {
        message.error("获取中药列表失败");
      }
    } catch (e: any) {
      message.error("获取中药列表失败," + e.message);
    }
  };

  useEffect(() => {
    loadData();
  }, [searchParams]);

  return (
    <PageContainer>
    <div>
      <Select defaultValue="chinese_name" onChange={value => setSearchType(value)}>
        <Option value="chinese_name">中文名称</Option>
        <Option value="english_name">英文名称</Option>
        <Option value="pinyin_name">拼音名称</Option>
      </Select>
      <Input
        placeholder="请输入搜索内容"
        onChange={e => setSearchValue(e.target.value.trim())}
        style={{ width: '200px', marginRight: '10px' }}
      />
      <Button
        type="primary"
        loading={loading}
        onClick={handleSearch}
      >
        搜索
      </Button>
    </div>
      <Table
        rowKey={"id"}
        pagination={{
          total: total,
          pageSizeOptions: [10, 20, 50, 100],
          onChange: loadData,
          showTotal: total => `共${total}条记录 `,
          defaultPageSize: 10,
          defaultCurrent: 1,
          position: ['bottomRight'],
          size: 'small',
          showSizeChanger: true,
          showQuickJumper: true,
        }}
        style={{marginTop: 20}}
        columns={columns}
        dataSource={tcmList}/>
    </PageContainer>
  );
};

export default Tcm;
