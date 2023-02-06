import multiprocessing
import time
from datetime import datetime
import os

class Process(multiprocessing.Process):
    def __init__(self, id):
        super(Process, self).__init__()
        self.id = id

class FiveSecProcess(Process):
    def run(self):
        print('Starting process {} at {}'.format(self.id, str(datetime.now())))
        time.sleep(5)
        print("Ending the process with id: {} at {}".format(self.id, str(datetime.now())))

class TenSecProcess(Process):
    def run(self):
        print('Starting process {} at {}'.format(self.id, str(datetime.now())))
        time.sleep(10)
        print("Ending the process with id: {} at {}".format(self.id, str(datetime.now())))

if __name__ == '__main__':
    print('No. of usable CPU cores: {}'.format(len(os.sched_getaffinity(0))))
    processes = []
    for i in range(0, 40, 2):
        processes.append(FiveSecProcess(i))
        processes.append(TenSecProcess(i + 1))
    for process in processes:
        process.start()
    for process in processes:
        process.join()