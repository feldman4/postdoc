#!/bin/bash
# to add completion, run with source: . app_completion.sh

app.sh -- --completion | sed "s/app.py/app.sh/g" > ~/.fire_completion
source ~/.fire_completion