#!/bin/bash

samtools view $1 $2 | awk '{sum+=$5} END { print "MAPQ=", sum/NR}'
