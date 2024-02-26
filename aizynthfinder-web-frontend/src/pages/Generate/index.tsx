import {PageContainer,} from '@ant-design/pro-components';
import React from 'react';
import {Button, Form, Input} from 'antd';

const onFinish = (values: any) => {
  console.log('Success:', values);
};

const onFinishFailed = (errorInfo: any) => {
  console.log('Failed:', errorInfo);
};

type FieldType = {
  smiles?: string;
};
const onChange = (e: React.ChangeEvent<HTMLInputElement | HTMLTextAreaElement>) => {
  console.log(e);
};
const TableList: React.FC = () => {


  return (
    <PageContainer>
      <Form
        name="basic"
        labelCol={{span: 8}}
        wrapperCol={{span: 16}}
        style={{maxWidth: 600}}
        initialValues={{remember: true}}
        onFinish={onFinish}
        onFinishFailed={onFinishFailed}
        autoComplete="off"
      >
        <Form.Item<FieldType>
          label="SMILES"
          name="smiles"
          rules={[{required: true, message: '该项不能为空!'}]}
        >
          <Input placeholder="请输入你想进行逆合成的SMILES式" allowClear onChange={onChange}/>
        </Form.Item>
        <Form.Item wrapperCol={{offset: 8, span: 16}}>
          <Button type="primary" htmlType="submit">
            确认
          </Button>
        </Form.Item>
      </Form>
    </PageContainer>
  );
};

export default TableList;
