#!/bin/bash
cd $(dirname $(readlink -f $0))
pytest parsers/hierarchicalconfigparser.py -s -v;
coverage run --source=. -m pytest parsers/hierarchicalconfigparser.py
coverage report -m parsers/hierarchicalconfigparser.py
#coverage report -m
