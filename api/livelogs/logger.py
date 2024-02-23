from datetime import datetime


class LiveLogger:
    def __init__(self):
        self.log_messages = []

    def push(self, message):
        self.log_messages.append(message)
        logging.info(message)


liveLogger = LiveLogger()