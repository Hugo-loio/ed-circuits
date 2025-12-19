#!/bin/sh

julia -e "using Pkg; Pkg.develop(path=\"$(dirname "$0")\")"
