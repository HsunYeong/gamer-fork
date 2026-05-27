#!/bin/bash

python3 construct_ic.py > >(tee -a ICs_Record__Note) 2> >(tee -a log_construct_ic >&2)
