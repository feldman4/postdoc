#!/bin/bash
# to add completion, run with source: . _appcompletion.sh

app.sh -- --completion \
 | sed "s/app.py/app.sh/g" \
 > ~/.fire_completion

app.sh fix_fire_completion ~/.fire_completion

source ~/.fire_completion
