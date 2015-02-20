#!/usr/bin/env python

import json
import sys
import math

def read_json_file(filename):
    with open(filename, "r") as file:
        return json.loads(file.read())

if __name__ == '__main__':
    if len(sys.argv) != 3:
        exit("Usage: %s <json> <field>" % sys.argv[0])
    else:
        data = read_json_file(sys.argv[1])
        print(data[sys.argv[2]])
