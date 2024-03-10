import io

import pandas as pd
import requests

from aizynthfinder.aizynthfinder import AiZynthFinder
from rdkit import Chem
from rdkit.Chem.Draw import MolToImage

from aizynthfinder.chem import Molecule
from aizynthfinder.context.scoring import CombinedScorer
from utils.HashUtils import HashUtils
from Config import Config


class GenerateService:
    @staticmethod
    def generate_img_from_smiles(smiles):
        res = None
        try:
            mol = Chem.MolFromSmiles(Molecule(smiles=smiles).smiles)
            img = MolToImage(mol)
            name = HashUtils.md5_hash(smiles)
            img.save(Config.FRONT_END_PATH + f"/public/images/{name}.png")
            res = name
        except:
            print(f"生成失败,无效的smiles表达式:[{smiles}]")
        return res

    @staticmethod
    def generate_route_from_smiles(smiles, scorers, finder):
        res = None
        try:
            # # 初始化 AiZynthFinder
            # filename = "../model_database/config.yml"
            # finder = AiZynthFinder(filename)
            # # 选择库存、扩展策略和过滤策略
            # finder.stock.select("zinc")
            # finder.expansion_policy.select("uspto_condition")
            # finder.filter_policy.select("uspto")
            # 设置目标 SMILES
            finder.target_smiles = smiles
            # 执行树搜索
            finder.tree_search()
            finder.build_routes()
            name = HashUtils.md5_hash(smiles)
            scorer = CombinedScorer(finder.config, scorers)
            print("评分函数初始化成功")
            trees = finder.routes.reaction_trees
            res_dict = {'routes': []}
            for i in range(len(trees)):
                trees[i].to_image().save(
                    Config.FRONT_END_PATH + f"/public/images/{name}{i + 1}.png")
                temp = trees[i].to_dict()
                temp['image_name'] = f"{name}{i + 1}"
                temp['index'] = i
                temp['score'] = round(scorer(trees[i]), 8)
                temp_dict = res_dict['routes']
                temp_dict.append(temp)
                res_dict['routes'] = temp_dict
            res_dict['routes'] = sorted(res_dict['routes'], key=lambda x: x['score'], reverse=True)
            res = res_dict
        except:
            print(f"逆合成失败:[{smiles}]")
        return res

    @staticmethod
    def get_mols_from_smiles(smiles):
        res = None
        try:
            url = "http://8.140.51.36:5000/mol2mol"
            data = {
                "smi_input": smiles
            }
            response = requests.post(url, json=data)
            result = response.json()
            csv_content = result['csv_content']
            string_io = io.StringIO(csv_content)
            data = pd.read_csv(string_io)
            length = min(len(data), 20)
            res_dict = {'mols': []}
            for i in range(length):
                temp = {}
                mol = Chem.MolFromSmiles(Molecule(smiles=data['SMILES'][i]).smiles)
                img = MolToImage(mol)
                name = HashUtils.md5_hash(data['SMILES'][i])
                img.save(Config.FRONT_END_PATH + f"/public/images/{name}.png")
                temp['key'] = i + 1
                temp['smiles'] = data['SMILES'][i]
                temp['tanimoto'] = data['Tanimoto'][i]
                temp['nll'] = data['NLL'][i]
                temp['image_name'] = name
                temp['generate_with_smiles'] = []
                temp_dict = res_dict['mols']
                temp_dict.append(temp)
                res_dict['mols'] = temp_dict
            res = res_dict
        except:
            print(f"分子生成失败:[{smiles}]")
        return res
