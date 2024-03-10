import {PageContainer,} from '@ant-design/pro-components';
import React, {useState} from 'react';
import {Button, Card, Divider, Form, Image, Input, message, Spin, Table} from 'antd';
import {getImgBySmiles, getMolsBySmiles, getRoutesBySmiles} from "@/services/aigenerate/GenerateController";

let dataSource: any[] = [];

const Generate: React.FC = () => {
  const [smilesToGen, setSmilesToGen] = useState<string>();
  const [originSmiles, setOriginSmiles] = useState<string>();
  const [submitting, setSubmitting] = useState<boolean>();
  const [loading, setLoading] = useState<boolean>();
  const [genAble, setGenAble] = useState<boolean>();
  const [expandedDataSource, setExpandedDataSource] = useState<{ [key: string]: any }>({});


  const selectToGenerate = async (smiles: string, key: number) => {
    console.log("当前操作的key为：" + key)
    if (key.toString() in expandedDataSource) {
      return;
    }
    const scorers: string[] = ['state score'];
    const param = {
      smiles,
      scorers,
    };
    let res = null;
    try {
      res = await getRoutesBySmiles(param, {});
      if (!res?.data) {
        message.error('无效的SMILES表达式:[' + smiles + '] -> [' + res.data + ']');
      } else {
        message.success('逆合成成功');
        const dict = expandedDataSource;
        dict[key.toString()] = res.data.routes;
        setExpandedDataSource(dict);
        console.log(res.data.routes)
        const mainTable = document.getElementById('main-table');
        if (mainTable) {
          console.log("获取到主表格")
          const mainTableComponent = (mainTable as any).reactInstance;
          if (mainTableComponent) {
            mainTableComponent.forceUpdate();
          }
        }
      }
      setLoading(false);
    } catch (e: any) {
      console.log(res)
      message.error('分析失败,' + e.message);
    }
  }

  const onFinish = async (values: any) => {
    dataSource = []
    setExpandedDataSource({})
    //避免重复提交
    if (submitting) {
      return;
    }
    setSubmitting(true);
    setOriginSmiles(values.smiles);
    const param = {
      ...values,
    };
    let res = null
    try {
      res = await getImgBySmiles(param, {});
      if (!res?.data) {
        message.error('无效的SMILES表达式:[' + values.smiles + '] -> [' + res.data + ']');
      } else {
        message.success('产物初始化成功,可以开始生成');
        setSmilesToGen('images/' + res.data + '.png');
        setGenAble(true);
        console.log(res.data)
      }
      setSubmitting(false);
    } catch (e: any) {
      console.log(res)
      message.error('分析失败,' + e.message);
      setSubmitting(false);
    }
  };
  const doGender = async () => {
    const values = {'smiles': originSmiles}
    console.log('尝试逆合成' + values.smiles)
    //避免重复提交
    if (loading) {
      return;
    }
    setLoading(true);
    const param = {
      ...values,
    };
    let res = null;
    try {
      res = await getMolsBySmiles(param, {});
      if (!res?.data) {
        message.error('无效的SMILES表达式:[' + values + '] -> [' + res.data + ']');
      } else {
        message.success('生成成功');
        console.log(res.data.mols)
        dataSource = res.data.mols;
      }
      setLoading(false);
    } catch (e: any) {
      console.log(res)
      message.error('生成失败,' + e.message);
      setLoading(false);
    }
  }

  const onFinishFailed = (errorInfo: any) => {
    console.log('Failed:', errorInfo);
  };
  const onSmilesChange = (e: React.ChangeEvent<HTMLInputElement | HTMLTextAreaElement>) => {
    console.log(e);
    setSmilesToGen('error');
    setGenAble(false);
    dataSource = []
    setExpandedDataSource({})
    setLoading(false);
  };
  const expandedColumns = [
    {
      title: '编号',
      dataIndex: 'index',
      key: 'index',
      align: 'center' as 'center',
      render: (text: any, record: any) => <>{record.index + 1}</>
    },
    {
      title: '合成路径',
      dataIndex: 'image_name',
      key: 'image_name',
      align: 'center' as 'center',
      render: (text: any, record: any) => <Image src={`images/${record.image_name}.png`} height={250}
                                                 alt={record.image_name}/>,
    },
    {
      title: '得分',
      dataIndex: 'score',
      key: 'score',
      align: 'center' as 'center',
    },
  ];
  const columns = [
    {
      title: '编号',
      dataIndex: 'key',
      key: 'key',
      align: 'center' as 'center',
    },
    {
      title: '分子smiles式',
      dataIndex: 'smiles',
      key: 'smiles',
      align: 'center' as 'center',
    },
    {
      title: '生成分子图像',
      dataIndex: 'image_name',
      key: 'image_name',
      align: 'center' as 'center',
      render: (text: any, record: any) => <Image src={`images/${record.image_name}.png`} width={350}
                                                 alt={record.image_name}/>,
    },
    {
      title: 'Tanimoto',
      dataIndex: 'tanimoto',
      key: 'tanimoto',
      align: 'center' as 'center',
    },
    {
      title: 'NLL',
      dataIndex: 'nll',
      key: 'nll',
      align: 'center' as 'center',
    },
  ];
  return (
    <PageContainer>
      <Card>
        <Form
          name="basic"
          initialValues={{remember: true}}
          onFinish={onFinish}
          onFinishFailed={onFinishFailed}
          autoComplete="off"
        >
          <Form.Item
            label="SMILES"
            name="smiles"
            rules={[{required: true, message: '该项不能为空!'}]}
          >
            <Input placeholder="请输入你想进行生成的SMILES式" allowClear onChange={onSmilesChange}/>
          </Form.Item>
          <Form.Item wrapperCol={{span: 16}}>
            <Button type="primary" htmlType="submit">
              确认传入分子
            </Button>
            <Divider type="vertical"/>
            <Button type="primary" disabled={!genAble} onClick={doGender} loading={loading}>
              开始生成
            </Button>
          </Form.Item>
          <Divider/>
          <Image
            title="传入分子"
            width={200}
            height={200}
            src={smilesToGen}
            fallback="data:image/png;base64,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"
          />
        </Form>
      </Card>
      {
        dataSource.length > 0 && (
          <Card title={`生成结果（共${dataSource.length}条）`}>
            <Table
              id="main-table"
              columns={columns}
              pagination={{pageSize: 10}}
              expandable={{
                onExpand: ((expanded, record) => {
                  selectToGenerate(record.smiles, record.key).then(() =>
                    setExpandedDataSource(expandedDataSource))
                }),
                expandedRowRender: (record) =>
                  <Card title="逆合成结果">
                    {
                      record.key.toString() in expandedDataSource && (
                        <Table
                          columns={expandedColumns}
                          dataSource={expandedDataSource[record.key.toString()]}
                          pagination={{pageSize: 10}}
                          expandable={{
                            showExpandColumn: false,
                          }}
                          scroll={{y: 800}}
                        />
                      )
                    }
                    {
                      !(record.key.toString() in expandedDataSource) && (
                        <Card>
                          <Spin/>
                        </Card>)
                    }
                  </Card>
              }}
              dataSource={dataSource}
            />
          </Card>
        )
      }
      {
        loading && (
          <Card title={`生成结果`}>
            <Spin/>
          </Card>
        )
      }
    </PageContainer>
  );
};

export default Generate;
