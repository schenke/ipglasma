#!/usr/bin/env python3
"""Doxygen input filter for the Markdown pages (README.md, OUTPUT.md, ...).

GitHub renders inline math written as $...$, Doxygen only \\f$...\\f$. This
filter rewrites $...$ into \\f$...\\f$ outside code spans and fenced code
blocks, so the same Markdown files render on GitHub and in the Doxygen
documentation. It follows GitHub's rule for inline math: the opening $ is
not preceded by a letter, digit or backslash and is followed by a
non-space, and the closing $ follows a non-space and is not followed by a
letter or digit.

Doxygen calls it with the file name (FILTER_PATTERNS in doc/Doxyfile.in)
and reads the filtered text from standard output.
"""

import re
import sys

INLINE_MATH = re.compile(r"(?<![\\$\w])\$(?=\S)([^$\n]*?\S)\$(?![\w$])")
CODE_SPAN = re.compile(r"(`+[^`]*`+)")


def convert_line(line):
    """Rewrites the inline math of one line, leaving its code spans."""
    parts = CODE_SPAN.split(line)
    for i in range(0, len(parts), 2):
        parts[i] = INLINE_MATH.sub(r"\\f$\1\\f$", parts[i])
    return "".join(parts)


def convert(text):
    """Rewrites the inline math of a Markdown text outside code blocks."""
    out = []
    in_fence = False
    for line in text.splitlines(keepends=True):
        if line.lstrip().startswith("```"):
            in_fence = not in_fence
            out.append(line)
        elif in_fence:
            out.append(line)
        else:
            out.append(convert_line(line))
    return "".join(out)


if __name__ == "__main__":
    with open(sys.argv[1], encoding="utf-8") as markdown:
        sys.stdout.write(convert(markdown.read()))
