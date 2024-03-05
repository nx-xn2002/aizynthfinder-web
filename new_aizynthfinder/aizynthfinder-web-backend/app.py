from datetime import timedelta

from flask import Flask, request, session
from flask_cors import CORS

from GenerateService import GenerateService
from models.BaseResponse import BaseResponse
from models.User import User

app = Flask(__name__)
app.secret_key = 'nx'
app.config['PERMANENT_SESSION_LIFETIME'] = timedelta(days=7)  # 设置为7天有效
app.config['SESSION_COOKIE_NAME'] = "aizynthfinder"
CORS(app, supports_credentials=True, origins='*', exposed_headers='Set-Cookie')


@app.route('/generate/getImgBySmiles', methods=['GET'])
def get_img_by_smiles():
    smiles = request.args.get('smiles')
    if smiles is not None:
        img = GenerateService.generate_img_from_smiles(smiles)
        return BaseResponse.success(img)
    else:
        return BaseResponse.null_error()


@app.route('/generate/getRoutesBySmiles', methods=['GET'])
def get_routes_by_smiles():
    scorers = request.args.getlist('scorers[]')
    print(scorers)
    smiles = request.args.get('smiles')
    if smiles is not None:
        res = GenerateService.generate_route_from_smiles(smiles, scorers)
        return BaseResponse.success(res)
    else:
        return BaseResponse.null_error()


@app.route('/user/login', methods=['POST'])
def login():
    username = request.json['username']
    password = request.json['password']
    if username is None or password is None:
        return BaseResponse.null_error()
    if username == "admin" and password == 'admin':
        user = User(username=username, avatar="logo.png")
        session['user'] = user
        print(f'用户{user}登录成功')
        return BaseResponse.success(user)
    else:
        return BaseResponse.params_error()


@app.route('/user/current', methods=['GET'])
def current():
    if 'user' in session:
        user = session.get('user')
        print(f'当前登录用户{user}')
        return BaseResponse.success(user)
    else:
        print('Not logged in')
        return BaseResponse.not_login()


@app.route('/user/logout', methods=['POST'])
def logout():
    user = session.pop('user', None)
    print(f'用户{user}退出登录')
    return BaseResponse.success("Logged out successfully!")


if __name__ == '__main__':
    app.run(host='localhost', port=9000, debug=True)
