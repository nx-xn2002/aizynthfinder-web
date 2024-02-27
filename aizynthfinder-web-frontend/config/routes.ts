export default [
  {
    path: '/user',
    layout: false,
    routes: [{name: '登录', path: '/user/login', component: './User/Login'}],
  },
  {path: '/welcome', name: '欢迎', icon: 'smile', component: './Welcome'},
  {name: 'SMILES逆合成', icon: 'DotChart', path: '/generate', component: './Generate'},
  {path: '/', redirect: '/welcome'},
  {path: '*', layout: false, component: './404'},
];
