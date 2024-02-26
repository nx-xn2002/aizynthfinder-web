from flask import Flask, request

from models.BaseResponse import BaseResponse

app = Flask(__name__)


@app.route('/')
def hello_world():
    response = BaseResponse.success("hello world")
    return response


@app.route('/generate/getImgBySmiles', methods=['GET'])
def get_img_by_smiles():
    smiles = request.args.get('smiles')
    print(smiles)
    if smiles is not None:
        return BaseResponse.success(smiles)
    else:
        return BaseResponse.null_error()


if __name__ == '__main__':
    app.run(host='0.0.0.0', port=10001, debug=True)
