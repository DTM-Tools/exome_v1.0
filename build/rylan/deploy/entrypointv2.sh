#!/bin/sh
sleep 2

# Run RyLAN
export RYLAN_DIR='/rylan'
cd /rylan
python3 $RYLAN_DIR/dtmtool/main_docker.py
