class BaseResponse(dict):
    def __init__(self, code=None, description=None, data=None):
        super(BaseResponse, self).__init__()
        self['code'] = code
        self['description'] = description
        self['data'] = data

    @staticmethod
    def success(data):
        return BaseResponse(0, "OK", data)

    @staticmethod
    def params_error():
        return BaseResponse(40000, "请求参数错误", None)

    @staticmethod
    def null_error():
        return BaseResponse(40000, "请求参数为空", None)

    @staticmethod
    def not_login():
        return BaseResponse(40100, "未登录", None)

    @staticmethod
    def no_auth():
        return BaseResponse(40101, "无权限", None)

    @staticmethod
    def system_error():
        return BaseResponse(50000, "系统内部异常", None)
