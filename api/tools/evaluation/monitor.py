
# pip install nvidia-ml-py
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
    

    def get_sys_info(self) -> dict:
        sys_info: dict = {}
        cpu, ram = self.get_cpu_mem_info()
        gpus = self.get_nvidia_info()
        sys_info['CPU'] = cpu
        sys_info['RAM'] = ram
        if len(gpus['gpus']) > 0:
            gpu_list = []
            for i in range(len(gpus['gpus'])):
                gpu = f"{gpus['gpus'][i]['gpu_model']} @ {gpus['gpus'][i]['total']} GB"
                gpu_list.append(gpu)
            
            sys_info['GPU'] = gpu_list

        return sys_info
    

    def get_cpu_mem_info(self):
        import platform

        n_cores = psutil.cpu_count(logical=False)
        # n_thread = psutil.cpu_count()
        freq = psutil.cpu_freq().current / 1000 # GHz
        cpu_model = platform.processor()
        mem_total = round(psutil.virtual_memory().total / 1024 / 1024 / 1024, 2) # GB
        # mem_free = round(psutil.virtual_memory().available / 1024 / 1024 / 1024, 2) # GB
        # mem_process_used = round(psutil.Process(os.getpid()).memory_info().rss / 1024 / 1024 / 1024, 2) # GB
        cpu = f"{cpu_model} {n_cores}-core @ {freq} GHz"
        ram = f"{mem_total} GB"

        return cpu, ram
    

    def get_nvidia_info(self):
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
                    "gpu_model": nvmlDeviceGetName(handle),
                    "total": round(memory_info.total / 1024 / 1024 / 1024, 2), # GB
                    "free": round(memory_info.free / 1024 / 1024 / 1024, 2), # GB
                    "used": round(memory_info.used / 1024 / 1024 / 1024, 2), # GB
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
        info = self.get_nvidia_info()
        if len(info['gpus']) > 0:
            used = info['gpus'][0]['used']
            tot = info['gpus'][0]['total']
            mem_rate = used/tot

        return mem_rate
