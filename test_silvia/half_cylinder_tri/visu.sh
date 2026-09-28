#!/bin/bash
set -e

ORDER=${1:-1}

cd outputs/o${ORDER}

pvbatch ../../visualize_fields.py

date > visu.timestamp
