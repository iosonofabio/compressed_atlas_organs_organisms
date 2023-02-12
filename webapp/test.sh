#!/bin/bash

# Adapted from python-igraph, original author Tamas Nepus:
# https://github.com/igraph/python-igraph/blob/709e7023aef4f4c4c93d385f4ed11adab6f7cbae/test.sh

###############################################################################

set -e

CLEAN=0
VENV_DIR=.venv
VERBOSE=0
MAIN_FLASK_FILE=app.py
CERT=0
export FLASK_DEBUG=1

if [ x$CLEAN = x1 ]; then
    rm -rf ${VENV_DIR}
fi

FLASK=${VENV_DIR}/bin/flask
PIP=${VENV_DIR}/bin/pip
PYTEST=${VENV_DIR}/bin/pytest

if [ ! -d ${VENV_DIR} ]; then
    python -m venv ${VENV_DIR}
fi
${VENV_DIR}/bin/pip install -r requirements.txt

if [ x$CERT = x1 ]; then
  if [ x$VERBOSE = x1 ]; then
    echo "${FLASK} run --cert=adhoc"
  fi
  ${FLASK} run --cert=adhoc
else
  ${FLASK} run
fi
