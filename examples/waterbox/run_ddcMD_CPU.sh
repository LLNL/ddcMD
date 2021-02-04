#!/bin/bash

# 1. link the restart file to running directory
ln -s snapshot.mem/restart
# 2. run executable
./ddcMD_CPU_kras_Debug
