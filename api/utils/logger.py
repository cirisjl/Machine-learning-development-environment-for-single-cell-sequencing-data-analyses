from loguru import logger as lg
import sys
import asyncio



class AppLogger:
    def __init__(self):
        self.app_logger = lg

    def set_logger(self, log_path=sys.stderr, user_id=None, rotation='1 MB', retention='1 days', filter_type=None, format="{time:YYYY-MM-DD HH:mm:ss} | {level} | {extra[user_id]} | {message}", level='DEBUG'):
        """
        :param log_path: log file path
        :param filter_type: filter
        :param level: [TRACE, DEBUG, INFO, SUCCESS, WARNING, ERROR, CRITICAL]
        :return:
        """
        if user_id:
            log_path = './' + user_id +'.log'
            
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



async def log_reader(user_id, n=5) -> list:
    """Log reader

    Args:
        n (int, optional): number of lines to read from file. Defaults to 5.

    Returns:
        list: List containing last n-lines in log file with html tags.
    """
    log_file = './' + user_id +'.log'
    log_lines = []
    with open(log_file, "r") as file:
        for line in file.readlines()[-n:]:
            if line.__contains__("ERROR"):
                log_lines.append(f'<span class="text-red-500">{line}</span><br/>')
            elif line.__contains__("WARNING"):
                log_lines.append(f'<span class="text-yellow-400">{line}</span><br/>')
            elif line.__contains__("SUCCESS"):
                log_lines.append(f'<span class="text-green-500">{line}</span><br/>')
            else:
                log_lines.append(f"{line}<br/>")
        return log_lines
    

user_format = "{level} | {message}"
logger = AppLogger()