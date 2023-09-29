#!/bin/bash
setup_dir=$(dirname "$0")
docker build -t jrboyd/tap $setup_dir/../docker