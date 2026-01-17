#!/bin/sh
set -euo # "strict" mode

uv pip install \
	meson \
	ninja \
	jsonschema \
	sphinx
