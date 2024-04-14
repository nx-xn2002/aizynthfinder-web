from flask_sqlalchemy import SQLAlchemy
from datetime import datetime

db = SQLAlchemy()


class UserModel(db.Model):
    __tablename__ = 'user'
    id = db.Column(db.Integer, primary_key=True, autoincrement=True, comment='用户id')
    username = db.Column(db.String(100), nullable=False, comment='用户名')
    password = db.Column(db.String(500), nullable=False, comment='密码')
    join_time = db.Column(db.DateTime, default=datetime.now, comment='加入时间')
    status = db.Column(db.Boolean, default=True, comment='是否启用')
    # ForeignKey 默认注册为普通用户
    role_id = db.Column(db.Integer, db.ForeignKey('role.id'), default=2, comment='用户角色')
    # Relationship
    roles = db.relationship('RoleModel', backref=db.backref('users', lazy='dynamic'))

    def to_dict(self):
        return {
            'id': self.id,
            'username': self.username,
            'createTime': self.join_time.strftime('%Y-%m-%d %H:%M:%S'),
            'status': self.status,
            'roles': self.roles.role_name,
        }


class RoleModel(db.Model):
    __tablename__ = 'role'
    id = db.Column(db.Integer, primary_key=True, autoincrement=True, comment='角色id')
    role_name = db.Column(db.String(100), nullable=False, comment='角色名称')
    role_desc = db.Column(db.String(100), nullable=False, comment='角色描述')


class TcmIngredientRelationModel(db.Model):
    __tablename__ = 'tcm_ingredient_relation'
    id = db.Column(db.Integer, primary_key=True, autoincrement=True, comment='id')
    tcm_id = db.Column(db.String(100), nullable=False, comment='中药id')
    ingredient_id = db.Column(db.String(100), nullable=False, comment='化学物质id')


class TcmModel(db.Model):
    __tablename__ = 'tcm'
    id = db.Column(db.Integer, primary_key=True, autoincrement=True, comment='中药id')
    chinese_name = db.Column(db.String(100), nullable=False, comment='中文名称')
    pinyin_name = db.Column(db.String(100), nullable=False, comment='拼音名称')
    english_name = db.Column(db.String(100), nullable=False, comment='英文名称')
    use_part = db.Column(db.String(100), nullable=True, comment='使用部位')
    property_flavor = db.Column(db.String(100), nullable=True, comment='性味')
    channel_tropism = db.Column(db.String(100), nullable=True, comment='归经')
    effect = db.Column(db.String(100), nullable=True, comment='功效')
    indication = db.Column(db.String(100), nullable=True, comment='主治')
    sdfile = db.Column(db.String(100), nullable=True, comment='SD文件')
    ref_source = db.Column(db.String(100), nullable=True, comment='来源')

    def to_dict(self):
        return {
            'id': self.id,
            'chineseName': self.chinese_name,
            'pinyinName': self.pinyin_name,
            'englishName': self.english_name,
            'usePart': self.use_part,
            'propertyFlavor': self.property_flavor,
            'channelTropism': self.channel_tropism,
            'effect': self.effect,
            'indication': self.indication,
            'sdfile': self.sdfile,
            'refSource': self.ref_source,
        }


class IngredientModel(db.Model):
    __tablename__ = 'ingredient'
    id = db.Column(db.Integer, primary_key=True, autoincrement=True, comment='id')
    name = db.Column(db.String(100), nullable=False, comment='名称')
    synonyms = db.Column(db.String(100), nullable=False, comment='同义词')
    smiles = db.Column(db.Text, nullable=False, comment='分子式')
    inchi = db.Column(db.Text, nullable=False, comment='化合物标识')
    inchikey = db.Column(db.Text, nullable=False, comment='哈希化合物标识')

    def to_dict(self):
        return {
            'id': self.id,
            'name': self.name,
            'synonyms': self.synonyms,
            'smiles': self.smiles,
            'inchi': self.inchi,
            'inchikey': self.inchikey,
        }
