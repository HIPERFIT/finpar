#!/usr/bin/env python

import json
import sys
import math

def read_json_file(filename):
    with open(filename, "r") as file:
        return json.loads(file.read())

def generate_platform_mk(data):
    if "gpu" in data:
        gpu=data["gpu"]
        print("HAVE_GPU = 1")
        print("GPU_DEVICE_ID = %d" % gpu["device_id"])
        print("GPU_WARP = %d" % gpu["warp_size"])
        print("GPU_LG_WARP = %d" % math.log(gpu["warp_size"],2))
        print("GPU_LOCAL_MEM = %d" % gpu["local_memory"])
        print("GPU_CONST_MEM = %d" % gpu["constant_memory"])
        print("GPU_REGISTER_MEM = %d" % gpu["register_memory"])
        print("GPU_DEVICE_MEM = %d" % (gpu["device_memory"]/(1024*1024)))
        print("GPU_LOCAL_MEM_PER_TH = %d" % gpu["local_memory_per_thread"])
        print("GPU_NUM_CORES = %d" % gpu["num_cores"])
    else:
        print("HAVE_GPU = 0")
    print("NCORES = %d" % data["num_cores"])
    print("MEMORY = %d" % data["main_memory"])
    print("export")

if __name__ == '__main__':
    if len(sys.argv) != 2:
        exit("Usage: %s <configuration-json>" % sys.argv[0])
    else:
        data = read_json_file(sys.argv[1])
        generate_platform_mk(data)

