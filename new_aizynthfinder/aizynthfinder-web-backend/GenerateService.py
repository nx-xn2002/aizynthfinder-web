from aizynthfinder.aizynthfinder import AiZynthFinder
from rdkit import Chem
from rdkit.Chem.Draw import MolToImage

from aizynthfinder.chem import Molecule
from utils.HashUtils import HashUtils


class GenerateService:
    @staticmethod
    def generate_img_from_smiles(smiles):
        res = None
        try:
            mol = Chem.MolFromSmiles(Molecule(smiles=smiles).smiles)
            img = MolToImage(mol)
            name = HashUtils.md5_hash(smiles)
            img.save(f"D:/MyWebProject/aizynthfinder-web/aizynthfinder-web-frontend/public/images/{name}.png")
            res = name
        except:
            print(f"生成失败,无效的smiles表达式:[{smiles}]")
        return res

    @staticmethod
    def generate_route_from_smiles(smiles):
        res = None
        try:
            # 初始化 AiZynthFinder
            filename = "../model_database/config.yml"
            finder = AiZynthFinder(filename)
            # 选择库存、扩展策略和过滤策略
            finder.stock.select("zinc")
            finder.expansion_policy.select("uspto_condition")
            finder.filter_policy.select("uspto")
            # 设置目标 SMILES
            finder.target_smiles = smiles
            # 执行树搜索
            finder.tree_search()
            finder.build_routes()
            name = HashUtils.md5_hash(smiles)
            trees = finder.routes.reaction_trees
            res_dict = {'routes': []}
            for i in range(len(trees)):
                trees[i].to_image().save(
                    f"D:/MyWebProject/aizynthfinder-web/aizynthfinder-web-frontend/public/images/{name}{i + 1}.png")
                temp = trees[i].to_dict()
                temp['image_name'] = f"{name}{i + 1}"
                temp['index'] = i
                temp_dict = res_dict['routes']
                temp_dict.append(temp)
                res_dict['routes'] = temp_dict
            res = res_dict
        except:
            print(f"逆合成失败:[{smiles}]")
        return res
