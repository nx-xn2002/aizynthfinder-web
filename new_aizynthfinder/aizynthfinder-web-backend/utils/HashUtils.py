import hashlib


class HashUtils:
    @staticmethod
    def md5_hash(input_string):
        hash_obj = hashlib.md5()
        hash_obj.update(input_string.encode('utf-8'))
        hashed_str = hash_obj.hexdigest()
        return hashed_str
