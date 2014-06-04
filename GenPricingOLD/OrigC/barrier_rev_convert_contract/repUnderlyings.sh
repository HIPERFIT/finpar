#!/bin/sh

sed -e "s/underlyings\[\([0-9a-zA-Z_]*\)\]\[\([0-9a-zA-Z_]*\)\]/underlyings(\1,\2)/g"  < $1 >$2
