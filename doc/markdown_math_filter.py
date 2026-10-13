#!/usr/bin/env python3
"""Doxygen input filter for the Markdown pages (README.md, OUTPUT.md, ...).

GitHub renders inline math written as $...$, Doxygen only \\f$...\\f$. This
filter rewrites $...$ into \\f$...\\f$ outside code spans and fenced code
blocks, so the same Markdown files render on GitHub and in the Doxygen
documentation. It follows GitHub's rule for inline math, so it converts
exactly the math GitHub renders: the opening $ starts the line or follows
a space or "(", and is followed by a non-space; the closing $ follows a
non-space and is not followed by a letter, digit or underscore.

Code spans and fences follow CommonMark: a code span closes with a run of
as many backticks as it opened with, and a fence (``` or ~~~, at least
three) closes with a line of at least as many of the same character.

Doxygen calls it with the file name (FILTER_PATTERNS in doc/Doxyfile.in)
and reads the filtered text from standard output.
"""

import re
import sys

INLINE_MATH = re.compile(r"(?<![^\s(])\$(?=\S)([^$\n]*?\S)\$(?![\w$])")
CODE_SPAN = re.compile(r"(?<!`)(`+)(?!`).+?(?<!`)\1(?!`)")
FENCE = re.compile(r"[ \t]*(`{3,}|~{3,})")


def convert_line(line):
    """Rewrites the inline math of one line, leaving its code spans."""
    out = []
    start = 0
    for span in CODE_SPAN.finditer(line):
        out.append(INLINE_MATH.sub(r"\\f$\1\\f$", line[start:span.start()]))
        out.append(span.group(0))
        start = span.end()
    out.append(INLINE_MATH.sub(r"\\f$\1\\f$", line[start:]))
    return "".join(out)


def convert(text):
    """Rewrites the inline math of a Markdown text outside code blocks."""
    out = []
    fence = None
    for line in text.splitlines(keepends=True):
        match = FENCE.match(line)
        if fence is None:
            if match:
                fence = match.group(1)
                out.append(line)
            else:
                out.append(convert_line(line))
        else:
            if (match and match.group(1)[0] == fence[0]
                    and len(match.group(1)) >= len(fence)
                    and not line[match.end():].strip()):
                fence = None
            out.append(line)
    return "".join(out)


if __name__ == "__main__":
    with open(sys.argv[1], encoding="utf-8") as markdown:
        sys.stdout.write(convert(markdown.read()))
