

**注意事项：开发时请勿将new_aizynthfinder包整个提交至git，提交前请运行git status进行检查**

**注意事项：开发时请勿将new_aizynthfinder包整个提交至git，提交前请运行git status进行检查**

**注意事项：开发时请勿将new_aizynthfinder包整个提交至git，提交前请运行git status进行检查**



# aizynthfinder-web（开发中）

基于aizynthfinder化学逆合成算法包制作的web网站的前后端代码仓库





## 前端 - aizynthfinder-web-frontend

### 技术选型

- HTML + CSS + JavaScript 三件套
- React
- Ant Design组件库
- Umi
- Ant Design Pro 脚手架



### 基本页面

#### 用户登录页

<img src="images/image-20240228040431631.png" height="400" style="float:left;">

#### 欢迎页

<img src="images/image-20240228040511349.png" height="400" style="float:left;">

#### 逆合成分析页

<img src="images/image-20240228185310364.png" height="400" style="float:left;">

## 后端 - aizynthfinder-web-backend

### 技术选型

- Python + Flask







## 快速开始

### 1. 将项目克隆至本地

```bash
git clone https://github.com/nx-xn2002/aizynthfinder-web.git
```
进入项目文件夹，创建名为`.git`的空文件夹，否则后面会报个错

### 2. 安装相关依赖

安装最新LTS版本的node.js以便运行前端项目 https://nodejs.org/en

在命令行更新当前npm至最新版本

### 3. 配置并启动前端项目

在idea、webstorm或pycharm等ide中，打开`aizynthfinder-web/aizynthfinder-web-frontend`文件夹，前端项目文件夹如图所示

<img src="images/image-20240302145112299.png" height="500" style="float:left;">

在当前`aizynthfinder-web-frontend`目录下，命令行运行`npm install`指令安装依赖

运行完毕后，运行`package.json`文件中的`start`命令启动前端（或者执行命令`npm run start`）

![image-20240302145719677](images/image-20240302145719677.png)

启动成功可在浏览器看到以下页面

<img src="images/image-20240302150829448.png" height="300" style="float:left;">

### 4. 拉取后端

在外部其他目录下，运行拉取后端依赖包

```bash
git clone https://github.com/chengfengke/new_aizynthfinder.git
```

然后在包目录下运行拉取大文件

```bash
git lfs pull
```

完成后，为避免项目中出现冲突，请将当前拉取到的`new_aizynthfinder`文件夹内的`.git`隐藏文件夹删除并重新创建一个名为`.git`的文件夹

严格完成以上步骤后，将当前文件夹拖入`aizynthfinder-web`文件夹中，然后在pycharm中打开`aizynthfinder-web/new_aizynthfinder`文件夹

### 5.后端环境配置

命令行运行指令初始化conda虚拟环境（需已安装anaconda或者其他工具）
```bash
conda create "python>=3.8,<3.10" -n aizynth-env
```

激活虚拟环境并安装相应依赖
```bash
conda activate aizynth-env
cd ./aizynthfinder-web-backend
pip install -r requirements.txt  #后续若加了新包记得同步更新到此txt文件中
```

在`aizynthfinder-web-backend/Config.py`文件中，修改`FRONT_END_PATH`的属性值，配置当前前端的路径，如：
```config
FRONT_END_PATH = r"E:\Projects\aizynthfinder-web\aizynthfinder-web-frontend"
```

修改`aizynthfinder-web/new_aizynthfinder/model_database/config_dev.yml`中的路径，可以使用`Ctrl+H`来批量替换为你本地真实的路径。如：
<img width="939" alt="32628f8ca8afe956aec6360e3712679" src="https://github.com/user-attachments/assets/2d7cb43c-5b7b-466b-b73e-fd60e1db1e81">

添加Python解释器，打开`Pycharm-File-Settings-Project:new_aizynthfinder-Python Interpreter-Add Interpreter-Add Local Interpreter-Conda Environmen-Use existing environment`，选择aizynth-env，点击OK，等待项目刷新。

配置好解释器以后Pycharm终端自动激活aizynth-env环境，不需要手动`conda activate`

### 6. 启动后端
接下来，启动`new_aizynthfinder/aizynthfinder-web-backend/app.py`即可

若被启动的是pytest in app，遵照如下解决方法：
![image](https://github.com/user-attachments/assets/1b749b71-48a5-4cdb-954f-d2255d61c7bc)

点击左上角加号，选择python

![image](https://github.com/user-attachments/assets/f3c01912-b229-4700-b841-bd435ec9356c)

路径按实际路径配置

![image](https://github.com/user-attachments/assets/1be32f24-c811-4f2d-961e-fee1fab54cdd)

启动新配置的app.py即可

在浏览器中使用默认账号密码admin尝试进行登录，若成功，则前后端配置启动完毕


## 开发
**前端中，配置前端请求后端地址的文件路径**：`\aizynthfinder-web\aizynthfinder-web-frontend\src\app.tsx`

**后端中，配置后端部署地址的文件路径**：`\aizynthfinder-web\new_aizynthfinder\aizynthfinder-web-backend\app.py`

二者应当保持一致

## 调试

**前端**：浏览器中按`F12`，可查看网页发出的请求与响应

**后端**：在pycharm的终端上可以看到输出的关于前端的请求记录
