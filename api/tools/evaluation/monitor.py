
import os
import time
import psutil
# import GPUtil
from threading import Thread
from pynvml import *



class Monitor(Thread):
    def __init__(self, delay=1):
        super(Monitor, self).__init__()
        self.stopped = False
        self.delay = delay # Time between calls to GPUtil
        self.time_points = []
        self.cpu_usage = []
        self.mem_usage = []
        self.gpu_mem_usage = []
        self.start()


    def run(self):
        while not self.stopped:
            # Obtaining all the essential details
            self.time_points.append(time.time())
            self.cpu_usage.append(psutil.cpu_percent())
            self.mem_usage.append(psutil.virtual_memory().percent)
            self.gpu_mem_usage.append(self.gpu_mem_percent())
            time.sleep(self.delay)
        

    def stop(self):
        self.stopped = True
        return self.time_points, self.cpu_usage, self.mem_usage, self.gpu_mem_usage
    

    def nvidia_info(self):
        # pip install nvidia-ml-py
        nvidia_dict = {
            "state": True,
            "nvidia_version": "",
            "nvidia_count": 0,
            "gpus": []
        }
        try:
            nvmlInit()
            nvidia_dict["nvidia_version"] = nvmlSystemGetDriverVersion()
            nvidia_dict["nvidia_count"] = nvmlDeviceGetCount()
            for i in range(nvidia_dict["nvidia_count"]):
                handle = nvmlDeviceGetHandleByIndex(i)
                memory_info = nvmlDeviceGetMemoryInfo(handle)
                gpu = {
                    "gpu_name": nvmlDeviceGetName(handle),
                    "total": memory_info.total,
                    "free": memory_info.free,
                    "used": memory_info.used,
                    "temperature": f"{nvmlDeviceGetTemperature(handle, 0)}℃",
                    "powerStatus": nvmlDeviceGetPowerState(handle)
                }
                nvidia_dict['gpus'].append(gpu)
        except NVMLError as _:
            nvidia_dict["state"] = False
        except Exception as _:
            nvidia_dict["state"] = False
        finally:
            try:
                nvmlShutdown()
            except:
                pass
        return nvidia_dict


    def gpu_mem_percent(self):
        mem_rate = 0.0
        info = self.nvidia_info()
        if len(info['gpus']) > 0:
            used = info['gpus'][0]['used']
            tot = info['gpus'][0]['total']
            mem_rate = used/tot

        return mem_rate
