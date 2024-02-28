import redis



class RedisClient:
    def __init__(self):
        # pool = redis.ConnectionPool(host='redis', port=6381, password='eYVX7EwVmmxKPCDmwMtyKVge8oLd2t81', db=0)
        # self.r = redis.Redis(connection_pool=pool)
        self.r = redis.Redis(host='redis', port=6381, password='eYVX7EwVmmxKPCDmwMtyKVge8oLd2t81', db=0)


    def write_log(self, user_id, msg):
        self.r.ltrim(user_id, 1, 0)
        self.r.lpush(user_id, msg)


    def read_log(self, user_id, start, end):
        log_lines = []
        logs = self.r.lrange(user_id, start, end)
        self.r.ltrim(user_id, start, end)

        for log in logs:
            log_lines.append(log.decode('utf-8'))

        return log_lines
    

    def clear_log(user_id, self):
        self.r.ltrim(user_id, 1, 0)


    def close(self):
        self.r.close()



redis_client = RedisClient()