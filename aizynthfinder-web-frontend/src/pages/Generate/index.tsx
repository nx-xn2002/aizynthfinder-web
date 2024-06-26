import {PageContainer,} from '@ant-design/pro-components';
import React, {useState} from 'react';
import {Button, Card, Divider, Form, Image, Input, message, Select, SelectProps, Spin, Table} from 'antd';
import {getImgBySmiles, getRoutesBySmiles} from "@/services/aigenerate/GenerateController";
import { useLocation } from 'umi';
import type { TableColumnsType, TableProps } from 'antd';

let dataSource: any[] = [];

const columns = [
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
// 分子式表格数据
type TableRowSelection<T> = TableProps<T>['rowSelection'];

interface DataType {
  key: React.Key;
  id: number;
  name: string;
  synonyms: string;
  smiles: string;
  inchi: string;
  inchikey: string;
}

const IngredientColumns: TableColumnsType<DataType> = [
  {
    title: '名称',
    dataIndex: 'name',
  },
  {
    title: '同义词',
    dataIndex: 'synonyms',
  },
  {
    title: '分子式',
    dataIndex: 'smiles',
  },
  {
    title: '化合物标识',
    dataIndex: 'inchi',
  },
  {
    title: '哈希化合物标识',
    dataIndex: 'inchikey',
  },
];


const Generate: React.FC = () => {
  const [originSmiles, setOriginSmiles] = useState<string>();
  const [form] = Form.useForm();
  const [smilesValue, setSmilesValue] = useState(''); // 定义状态变量
  const location = useLocation();
  //console.log(location.state);
  const data: DataType[] = (location.state && location.state.ingredients) || [];
  // 使用 map 函数遍历 ingredients 数组，并为每个对象添加 key 属性
  const dataWithKeys = data.map((ingredient, i) => ({
    ...ingredient, // 保留原始对象的属性
    key: i, // 添加 key 属性，赋值为当前索引 i
  }));
  //分子式选择
  const [selectedRowKeys, setSelectedRowKeys] = useState<React.Key[]>([]);
  const onSelectChange = (newSelectedRowKeys: React.Key[]) => {
    console.log('selectedRowKeys changed: ', newSelectedRowKeys);
    setSelectedRowKeys(newSelectedRowKeys);
    form.setFieldsValue({ smiles: dataWithKeys[newSelectedRowKeys[0]].smiles });
  };

  const rowSelection: TableRowSelection<DataType> = {
    type: 'radio',
    selectedRowKeys,
    onChange: onSelectChange,
  };

  const [smilesToGen, setSmilesToGen] = useState<string>();
  const [scorers, setScorers] = useState<string[]>();
  const [submitting, setSubmitting] = useState<boolean>();
  const [loading, setLoading] = useState<boolean>();
  const [genAble, setGenAble] = useState<boolean>();

  const options: SelectProps['options'] = [
    {label: "state score", value: "state score"},
    {label: "max transform", value: "max transform"},
    {label: "fraction in stock", value: "fraction in stock"},
    {label: "number of reactions", value: "number of reactions"},
    {label: "number of pre-cursors", value: "number of pre-cursors"},
    {label: "number of pre-cursors in stock", value: "number of pre-cursors in stock"},
    {label: "average template occurrence", value: "average template occurrence"},
    {label: "sum of prices", value: "sum of prices"},
    {label: "route cost", value: "route cost"},
    {label: "number of reactions", value: "number of reactions"},
  ];
  const handleChange = (value: string[]) => {
    setScorers(value);
    if (genAble && value.length <= 0) {
      setGenAble(false);
    } else if (smilesToGen !== "error") {
      setGenAble(true);
    }
  };

  const onFinish = async (values: any) => {
    dataSource = []
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
        message.success('产物初始化成功,可以开始逆合成');
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
    const values = {'smiles': originSmiles, 'scorers': scorers}
    console.log('尝试逆合成' + values.smiles + values.scorers)
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
      res = await getRoutesBySmiles(param, {});
      if (!res?.data) {
        message.error('无效的SMILES表达式:[' + values + '] -> [' + res.data + ']');
      } else {
        message.success('逆合成成功');
        console.log(res.data.routes)
        dataSource = res.data.routes;
      }
      setLoading(false);
    } catch (e: any) {
      console.log(res)
      message.error('分析失败,' + e.message);
      setLoading(false);
    }
  }

  const onFinishFailed = (errorInfo: any) => {
    console.log('Failed:', errorInfo);
  };
  const onSmilesChange = (e: React.ChangeEvent<HTMLInputElement | HTMLTextAreaElement>) => {
    setSmilesValue(e.target.value);
    console.log(e);
    setSmilesToGen('error');
    setGenAble(false);
    dataSource = []
    setLoading(false);
  };
  return (
    <PageContainer>
      <Card>
        <Form
          name="basic"
          form={form}
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
            <Input placeholder="请输入你想进行逆合成的SMILES式" allowClear onChange={onSmilesChange} value={smilesValue}/>
          </Form.Item>
          <Form.Item wrapperCol={{span: 16}}>
            <Button type="primary" htmlType="submit">
              确认目标产物
            </Button>
            <Divider type="vertical"/>
            <Button type="primary" disabled={!genAble} onClick={doGender} loading={loading}>
              开始逆合成
            </Button>
          </Form.Item>
          <Form.Item
            label="评分函数"
            name="scorers"
          >
            <Select
              mode="multiple"
              allowClear
              style={{width: '100%'}}
              placeholder="请选择至少一项评分函数"
              onChange={handleChange}
              options={options}
            />
          </Form.Item>
           <Table rowSelection={rowSelection} columns={IngredientColumns} dataSource={dataWithKeys} />
          <Divider/>
          <Image
            title="目标产物"
            width={200}
            height={200}
            src={smilesToGen}
            fallback="data:image/png;base64,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"
          />
        </Form>
      </Card>
      {
        dataSource.length > 0 && (
          <Card title={`逆合成结果（共${dataSource.length}条）`}>
            <Table
              columns={columns}
              dataSource={dataSource}
              pagination={{pageSize: 10}}
              expandable={{
                showExpandColumn: false,
              }}
            />
          </Card>
        )
      }
      {
        loading && (
          <Card title={`逆合成结果`}>
            <Spin/>
          </Card>
        )
      }
    </PageContainer>
  );
};

export default Generate;
