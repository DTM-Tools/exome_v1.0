# DTM-Tools exome v1.0:Software Description

## Introduction
DTM-Tools was designed for prediction of blood antigen phenotypes from large NGS datasets. It is written in Python and designed to work in a Linux/Unix environment. It can only be run from the command line or through a Docker container – graphical and web interfaces will be implemented in the future.

The description and tutorial material presented in this document is intended as a complement, not a repetition, for the software’s publication. Please refer to ```Citation pending``` before reviewing this material.

DTM-Tools is written in object-oriented notation as a series of stages (static methods) that create objects, execute, and return the output to the pipeline state. Logging levels are incorporated for future web server deployment.

## Source Code and Structure

[Github main]: https://github.com/DTM-Tools/exome_v1.0

DTM-Tools source code can be accessed through [Github][Github main]: https://github.com/DTM-Tools/exome_v1.0. Please refer to the DTM-Tools publication for further details and laboratory validation data.

The following four subdirectories are included in exome v1.0:

_/build_: contains Docker definition files and the Docker deployment entrypoint.

_/databases_: contains three .csv files that define genomic coordinates and their interpretation rules. Described in further detail below.

_/rylantool_: contains DTM-Tools app python scripts

_/queryTools_: contains sample pymongo query scripts


