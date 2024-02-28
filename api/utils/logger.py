from loguru import logger as lg
import sys
import os



class AppLogger:
    def __init__(self):
        self.app_logger = lg

    def set_logger(self, log_path=sys.stderr, user_id=None, rotation='500 MB', retention='7 days', filter_type=None, format="{time:YYYY-MM-DD HH:mm:ss} | {level} | {extra[user_id]} | {message}", level='DEBUG'):
        """
        :param log_path: log file path
        :param filter_type: filter
        :param level: [TRACE, DEBUG, INFO, SUCCESS, WARNING, ERROR, CRITICAL]
        :return:
        """
        if user_id:
            log_path = './' + user_id +'.log'

        if not os.path.exists(os.path.dirname(log_path)):
             os.makedirs(os.path.dirname(log_path))
            
        dic = dict(
            sink=log_path,
            rotation=rotation,
            retention=retention,
            format=format,
            encoding='utf-8',
            level=level,
            enqueue=True,
        )
        if user_id:
            dic["filter"] = lambda record: record["extra"].get("user_id") == user_id
            # self.app_logger = self.app_logger.bind(user_id=user_id)
        elif filter_type:
            dic["filter"] = lambda x: filter_type in str(x['level']).upper()
        
        self.app_logger.add(**dic)

        return self.app_logger

    @property
    def get_logger(self):
        return self.app_logger
    
    @staticmethod
    def trace(self, user_id, msg):
        self.app_logger.bind(user_id=user_id).trace(msg)
 
    def debug(self, user_id, msg):
        self.app_logger.bind(user_id=user_id).debug(msg)
 
    def info(self, user_id, msg):
        self.app_logger.bind(user_id=user_id).info(msg)
 
    def success(self, user_id, msg):
        self.app_logger.bind(user_id=user_id).success(msg)
 
    def warning(self, user_id, msg):
        self.app_logger.bind(user_id=user_id).warning(msg)
 
    def error(self, user_id, msg):
        self.app_logger.bind(user_id=user_id).error(msg)
 
    def critical(self, user_id, msg):
        self.app_logger.bind(user_id=user_id).critical(msg)


# user_format = "{level} | {message}"
logger = AppLogger()
logger.set_logger('/data/logs/error.log', filter_type='ERROR')
logger.set_logger('/data/logs/activity.log', filter_type='INFO', level='INFO')