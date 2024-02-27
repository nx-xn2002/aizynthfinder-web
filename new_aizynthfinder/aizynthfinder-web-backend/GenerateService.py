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
