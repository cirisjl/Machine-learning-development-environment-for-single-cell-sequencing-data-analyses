import redis
import asyncio
from utils.logger import *


# unique_id: user_id or task_id
class RedisLogger:
    def __init__(self):
        pool = redis.ConnectionPool(host='redis', port=6381, password='eYVX7EwVmmxKPCDmwMtyKVge8oLd2t81', db=0)
        self.r = redis.Redis(connection_pool=pool)
        # self.r = redis.Redis(host='redis', port=6381, password='eYVX7EwVmmxKPCDmwMtyKVge8oLd2t81', db=0)
        # self.r = redis.Redis(host='127.0.0.1', port=6379, db=0)
        self.app_logger = logger
        

    def write_log(self, unique_id, msg):
        self.r.lpush(unique_id, msg)


    def read_log(self, unique_id, start, end):
        log_lines = []
        logs = self.r.lrange(unique_id, start, end)
        self.r.ltrim(unique_id, start, end)

        for log in logs:
            log_lines.append(log.decode('utf-8'))

        self.r.expire(unique_id, 60*30) # key expires in 30 minutes

        return log_lines[::-1]
    

    def clear_log(self, unique_id):
        self.r.ltrim(unique_id, 1, 0)


    def close(self):
        self.r.close()

    
    @staticmethod
    def trace(self, unique_id, msg):
        self.write_log(unique_id, 'TRACE    | ' + msg)
        self.app_logger.trace(unique_id, msg)
 

    def debug(self, unique_id, msg):
        self.write_log(unique_id, 'DEBUG    | ' + msg)
        self.app_logger.debug(unique_id, msg)
 

    def info(self, unique_id, msg):
        self.write_log(unique_id, 'INFO     | ' + msg)
        self.app_logger.info(unique_id, msg)
 

    def success(self, unique_id, msg):
        self.write_log(unique_id, 'SUCCESS  | ' + msg)
        self.app_logger.success(unique_id, msg)
 

    def warning(self, unique_id, msg):
        self.write_log(unique_id, 'WARNING  | ' + msg)
        self.app_logger.warning(unique_id, msg)
 

    def error(self, unique_id, msg):
        self.write_log(unique_id, 'ERROR    | ' + msg)
        self.app_logger.error(unique_id, msg)
 

    def critical(self, unique_id, msg):
        self.write_log(unique_id, 'CRITICAL | ' + msg)
        self.app_logger.critical(unique_id, msg)



redislogger = RedisLogger()

async def log_reader(unique_id, start=0, end=30) -> list:
    log_lines = []
    for line in redislogger.read_log(unique_id, start, end):
        if line.__contains__("ERROR"):
            log_lines.append(f'<span class="text-red-500">{line}</span><br/>')
        elif line.__contains__("WARNING"):
            log_lines.append(f'<span class="text-yellow-400">{line}</span><br/>')
        elif line.__contains__("SUCCESS"):
            log_lines.append(f'<span class="text-green-500">{line}</span><br/>')
        else:
            log_lines.append(f"{line}<br/>")
    return log_lines