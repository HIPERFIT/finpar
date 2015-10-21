#!/usr/bin/env python

import json
import sys

def read_json_file(filename):
    with open(filename, "r") as file:
        return json.loads(file.read())

if __name__ == '__main__':
    if len(sys.argv) < 3:
        exit("Usage: %s <file> [fields...]" % sys.argv[0])
    else:
        data = read_json_file(sys.argv[1])
        for field in sys.argv[2:]:
            try:
                print(data[field])
            except KeyError:
                exit("Key %s not found in JSON file." % field)
