#!/bin/bash
setup_dir=$(dirname "$0")
docker build -t tap $setup_dir/../docker