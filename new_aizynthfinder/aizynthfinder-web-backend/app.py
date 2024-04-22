from datetime import timedelta

import sqlalchemy
import Config
from aizynthfinder.aizynthfinder import AiZynthFinder
from flask import Flask, request, session, jsonify
from flask_cors import CORS

from GenerateService import GenerateService
from models.BaseResponse import BaseResponse
from models.User import User
from database_models import *
from utils.chat_utils import *
from selfies import encoder, decoder
import argparse

app = Flask(__name__)
app.secret_key = 'nx'
app.config['PERMANENT_SESSION_LIFETIME'] = timedelta(days=7)  # 设置为7天有效
app.config['SESSION_COOKIE_NAME'] = "aizynthfinder"
CORS(app, supports_credentials=True, origins='*', exposed_headers='Set-Cookie')
app.config.from_object(Config)
db.init_app(app)


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
        res = GenerateService.generate_route_from_smiles(smiles, scorers, finder)
        return BaseResponse.success(res)
    else:
        return BaseResponse.null_error()


@app.route('/getMolsBySmiles', methods=['GET'])
def get_mols_by_smiles():
    smiles = request.args.get('smiles')
    if smiles is not None:
        res = GenerateService.get_mols_from_smiles(smiles)
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


@app.route('/tcm/search', methods=['POST'])
def get_tcms():
    # page = int(request.json['current'])  # 获取页码，默认为第一页
    page = int(request.json.get('current', 1))
    # 获取每页显示的数据量，默认为 10 条
    per_page = int(request.json.get('pageSize', 10))
    # 获取中文名，默认为空字符串
    chinese_name = request.json.get('chinese_name', '')
    # 获取英文名，默认为空字符串
    english_name = request.json.get('english_name', '')
    # 获取拼音名，默认为空字符串
    pinyin_name = request.json.get('pinyin_name', '')

    # 构造查询语句
    query = TcmModel.query  # 使用 TcmModel 模型进行查询
    if chinese_name:  # 如果name不为空
        query = query.filter(TcmModel.chinese_name.ilike(f'%{chinese_name}%'))  # 使用 ilike() 方法查询所有包含 username 的用户名
    if english_name:  # 如果name不为空
        query = query.filter(TcmModel.english_name.ilike(f'%{english_name}%'))
    if pinyin_name:  # 如果name不为空
        query = query.filter(TcmModel.pinyin_name.ilike(f'%{pinyin_name}%'))

    # 分页查询
    pagination = query.paginate(page=page, per_page=per_page, error_out=False)  # 使用 paginate() 方法进行分页查询，不抛出异常
    tcms = pagination.items  # 获取当前页的数据
    total = pagination.total  # 获取总数据量
    # 构造返回数据
    data = {
        'list': [tcm.to_dict() for tcm in tcms],  # 将当前页的所有数据转换为字典形式，并存储在列表中
        'total': total,  # 总数据量
    }
    return BaseResponse.success(data)


@app.route('/ingredient/search', methods=['POST'])
def get_ingredient():
    tcm_id = request.json.get('id', '')  # 获取 tcm_id

    # 构造子查询
    subquery = db.session.query(TcmIngredientRelationModel.ingredient_id). \
        join(TcmModel, TcmModel.id == TcmIngredientRelationModel.tcm_id). \
        filter(TcmModel.id == tcm_id).subquery()

    # 构造主查询
    query = db.session.query(IngredientModel). \
        filter(IngredientModel.id.in_(subquery))

    # 如果 tcm_id 存在，则添加过滤条件
    if tcm_id:
        query = query.filter(TcmModel.id == tcm_id)

    ingredients = query.all()  # 执行查询并获取所有结果

    data = {
        'list': [ingredient.to_dict() for ingredient in ingredients],
        'total': len(ingredients),  # 总数为结果列表的长度
    }
    return BaseResponse.success(data)


@app.route('/associate/evaluate', methods=['POST'])
def evaluate_route():
    # 定义命令行参数
    parser = argparse.ArgumentParser(description='Description of your program')
    parser.add_argument('--CLI', action='store_true', help='Description of CLI parameter')
    parser.add_argument('--protein', action='store_true', help='Description of protein parameter')
    parser.add_argument('--load_8bit', action='store_true', help='Description of load_8bit parameter')
    parser.add_argument('--base_model', type=str, help='Description of base_model parameter')
    parser.add_argument('--share_gradio', action='store_true', help='Description of share_gradio parameter')
    parser.add_argument('--lora_weights', type=str, help='Description of lora_weights parameter')
    # 设置默认值
    parser.set_defaults(CLI=True, protein=False, load_8bit=False, base_model='/root/autodl-tmp/Llama-2-7b-chat-hf',
                        share_gradio=True,
                        lora_weights='/root/autodl-tmp/llama2-molinst-molecule-7b')

    args = parser.parse_args()

    # 根据命令行参数实例化函数类
    generate = main(args.CLI, args.protein, args.load_8bit, args.base_model, args.share_gradio,
                    args.lora_weights)

    data = request.json
    instruction = data.get('instruction')
    input_text = data.get('input_text')
    temperature = data.get('temperature', 0.1)
    top_p = data.get('top_p', 0.75)
    top_k = data.get('top_k', 40)
    num_beams = data.get('num_beams', 4)
    repetition_penalty = data.get('repetition_penalty', 1)
    max_new_tokens = data.get('max_new_tokens', 128)
    print(instruction, input_text, temperature, top_p, top_k, num_beams, repetition_penalty, max_new_tokens)
    #
    if generate.is_smiles(input_text):
        input_text = encoder(input_text)

    output = generate.evaluate_instruction(instruction, input_text, temperature, top_p, top_k, num_beams,
                                           repetition_penalty,
                                           max_new_tokens)

    if generate.is_selfies(output):
        try:
            smiles_output = Chem.MolToSmiles(Chem.MolFromSmiles(decoder(output)))
            output = smiles_output if smiles_output else output
        except Exception as e:
            pass

    data = {
        'result': output,
    }
    return BaseResponse.success(data)


def test_database_connection():
    with app.app_context():
        with db.engine.connect() as conn:
            res = conn.execute(sqlalchemy.text('select 1'))
            if res.fetchone()[0] == 1:
                print('Database connection successful')
            else:
                print('Database connection failed')


if __name__ == '__main__':
    # 初始化 AiZynthFinder
    filename = "../model_database/config_dev.yml"
    finder = AiZynthFinder(filename)
    # 选择库存、扩展策略和过滤策略
    finder.stock.select("zinc")
    finder.expansion_policy.select("uspto_condition")
    finder.filter_policy.select("uspto")
    test_database_connection()
    app.run(host='127.0.0.1', port=9000, debug=False)
