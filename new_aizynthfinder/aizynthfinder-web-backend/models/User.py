class User(dict):
    def __init__(self, username=None, password=None, avatar=None):
        super(User, self).__init__()
        self['username'] = username
        self['password'] = password
        self['avatar'] = avatar
