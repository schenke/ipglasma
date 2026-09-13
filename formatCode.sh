#!/usr/bin/env bash

if command -v clang-format &> /dev/null
then
    echo "formatting code ..."
    find src tests utilities \( -iname '*.h' -o -iname '*.cpp' \) \
        -not -path '*/third_party/*' | xargs clang-format -i -style=file
else
    echo "clang-format not found, skip formatting code"
fi
