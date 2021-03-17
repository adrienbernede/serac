#!/bin/bash

if [[ ! ${1} ]]
then
    echo "ERROR: expecting reference for uberenv repo" >&2
else
    uberenv_ref="${1}"
fi

cd scripts/uberenv
git fetch origin
git checkout ${uberenv_ref}
