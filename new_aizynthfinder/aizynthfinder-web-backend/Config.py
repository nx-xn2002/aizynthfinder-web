class Config:
    FRONT_END_PATH = "/www/wwwroot/aizynthfinder-web/aizynthfinder-web-frontend"


# 数据库配置
# HOSTNAME = '127.0.0.1'
HOSTNAME = '39.105.129.140'
PORT = 12308
USERNAME = 'aidrug'
PASSWORD = 'AiDrug!123456'
DATABASE = 'ai_drug_research'
DB_URI = 'mysql+pymysql://{}:{}@{}:{}/{}?charset=utf8mb4'.format(USERNAME, PASSWORD, HOSTNAME, PORT, DATABASE)
# DB_URI = f'mysql+pymysql://{USERNAME}:{PASSWORD}@{HOSTNAME}:{PORT}/{DATABASE}?charset=utf8mb4'
SQLALCHEMY_DATABASE_URI = DB_URI
