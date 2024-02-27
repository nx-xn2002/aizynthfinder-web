from flask import Flask, request
from flask_cors import CORS

from GenerateService import GenerateService
from models.BaseResponse import BaseResponse
from models.User import User

app = Flask(__name__)
CORS(app, supports_credentials=True)


@app.route('/')
def hello_world():
    response = BaseResponse.success("hello world")
    return response


@app.route('/generate/getImgBySmiles', methods=['GET'])
def get_img_by_smiles():
    smiles = request.args.get('smiles')
    if smiles is not None:
        img = GenerateService.generate_img_from_smiles(smiles)
        return BaseResponse.success(img)
    else:
        return BaseResponse.null_error()


@app.route('/user/login', methods=['POST'])
def login():
    username = request.form.get('username')
    password = request.form.get('password')
    if username is None or password is None:
        return BaseResponse.null_error()
    if username == "admin" and password == 'admin':
        return BaseResponse.success(User(username="admin", avatar="logo.png"))
    else:
        return BaseResponse.params_error()


if __name__ == '__main__':
    app.run(host='0.0.0.0', port=8000, debug=True)
