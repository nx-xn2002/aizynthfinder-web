export default [
  {
    path: '/user',
    layout: false,
    routes: [{name: '登录', path: '/user/login', component: './User/Login'}],
  },
  {path: '/welcome', name: '欢迎', icon: 'smile', component: './Welcome'},
  {name: 'SMILES逆合成', icon: 'DotChart', path: '/generate', component: './Generate'},
  {name: '分子生成器', icon: 'UnorderedList', path: '/mol2mol', component: './MolToMol'},
  {name: '中药查询', icon: 'SearchOutlined', path: '/tcm', component: './Tcm'},
  {name: '助理咨询', icon: 'RedditOutlined', path: '/associate', component: './Associate'},
  {name: '分子属性预测', icon: 'DotChart', path: '/predict', component: './Predict'},
  {path: '/', redirect: '/welcome'},
  {path: '*', layout: false, component: './404'},
];
