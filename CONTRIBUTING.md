# Contributing

## Adding an input parameter

1. Add a documented field with an initializer to the matching group struct in
   `src/Parameters.h` (e.g. `JimwlkParameters`). Name it like its input-file
   key.
2. Add one entry to the table in `src/ParameterTable.cpp`, e.g.
   `param("size", &P::lattice, &LatticeParameters::size).check(even())`, with
   `.optional("<default>")` if the key may be omitted, `.onlyIf(...)` if it is
   only read in some configurations, and `.check(...)` for its valid range.
   The parameters are read in the order of the table, so an `.onlyIf()`
   condition can only use parameters listed before it. Checks that combine
   several parameters go into `Parameters::validationErrors()`.
3. Document it in the input-parameter list of the README, with "read with …"
   if it has an `.onlyIf()` condition, and add a CHANGELOG entry.

Reading, the unknown-key check and the `usedParameters` output all follow
from the table; nothing else in the code needs to change. When you rename or
remove a key, add the old key to `renamedKeys()` or `replacedKeys()` in
`src/ParameterTable.cpp`, so that an old input file gets a hint instead of
only "unknown parameter".

## Changing an output file

Every file IP-Glasma writes is described in [OUTPUT.md](OUTPUT.md): when it is
written, its name, layout, column order and units. When you add an output
file or change the content or layout of one, update its section there in the
same change, and point to it from the Doxygen comment of the function that
writes it.

## License headers and authors

IP-Glasma is licensed under the GNU General Public License, version 3 or
later (see `LICENSE`). Every C++, Python and shell file of the project starts
with the header (with `#` instead of `//` in scripts, after the `#!` line)

```cpp
// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2011-2026 The IP-Glasma authors (see AUTHORS)
// This file is part of IP-Glasma; the license text is in LICENSE.
```

Start a new file with this header, and add yourself to `AUTHORS` with your
first contribution. Code taken from another project keeps its own copyright
and license notice; list it at the end of `AUTHORS`.

## Code formatting

The code in `src/`, `tests/` and `utilities/` is formatted with clang-format
19, the version the GitHub Actions workflow uses (e.g.
`pip install clang-format==19.1.7`; older versions cannot read
`.clang-format`). To format it, run from the repository root:

```
find src tests utilities \( -iname '*.h' -o -iname '*.cpp' \) -not -path '*/third_party/*' | xargs clang-format -i -style=file
```

For a pull request from a branch of this repository,
`.github/workflows/clang-format.yml` runs the same command and pushes a
`style: apply clang-format` commit if anything changed; pull that commit
before you continue working on the branch. Pull requests from forks are not
formatted automatically, so format them before opening the pull request.

## Code documentation

IP-Glasma uses [Doxygen](https://www.doxygen.nl/)-style documentation,
following the conventions used by
[SMASH](https://github.com/smash-transport/smash) (see its
[CONTRIBUTING.md](https://github.com/smash-transport/smash/blob/main/CONTRIBUTING.md)).
Every class, function, and member variable in `src/` should be documented
this way -- see "Scope" below.

### Building the documentation

You need Doxygen installed (and, optionally, Graphviz's `dot` for the
class/collaboration/call graphs). From your build directory:

```
cmake --build . --target doc
```

and open `doc/html/index.html` in a browser. Two more targets check
completeness (both build the documentation with only fully-documented
entities extracted, so any gap shows up as a Doxygen "not documented"
warning):

```
cmake --build . --target undocumented        # lists every warning
cmake --build . --target undocumented_count  # just counts them
```

`.github/workflows/doxygen.yml` runs `undocumented_test` on every push/PR
to `devel`/`main`, so a newly-added undocumented class/function/member
fails CI; run it locally first (see above) to catch this before pushing.
There is currently no hosted/published copy of the generated documentation
itself -- that's a separate, later step.

### Comment style

Use a `/** ... */` block, with continuation lines starting with a
leading `*`:

```cpp
/**
 * One-sentence summary, used as the brief description.
 *
 * Optional longer explanation, e.g. the physics/formula behind the
 * quantity, or which stage of the pipeline sets/consumes it.
 * \param[in] x What x is [unit].
 * \return What is returned [unit].
 */
```

- Do **not** add an explicit `\brief` tag for a plain one-sentence summary:
  the first sentence (up to the first period) is automatically used as the
  brief description. Use `\brief` only if you need the summary to end
  before a period that isn't the end of the sentence (rare).
- Blank comment line between the summary and a longer explanation, same as
  a blank line between paragraphs in prose.

### Tags

- `\param[in]`, `\param[out]`, `\param[in,out]` -- one per parameter, with
  the direction that actually applies. Every parameter gets one.
- `\return` -- required for every non-`void` return.
- `\note` / `\warning` -- for a real caveat (a non-obvious precondition, a
  surprising edge case), not for restating what the code already says.
- `\tparam` -- for template parameters.

### Units and physics notation

- Every documented physical quantity states its unit in trailing square
  brackets: `[GeV]`, `[fm]`, `[lattice units]`, `[dimensionless]`. This
  codebase mixes all of these freely (often for the same quantity at
  different stages), so never leave it implicit.
- Inline math uses `\f$ ... \f$` (e.g. `\f$T^{\tau\tau}\f$`,
  `\f$u^\mu\f$`); a displayed equation uses `\f[ ... \f]`.
- Reference the Milne-coordinate convention this codebase uses
  (\f$\tau, x, y, \eta\f$) explicitly wherever it disambiguates a
  quantity (e.g. an energy-momentum tensor component).

### Publication references

- Cite papers with `\cite <INSPIRE texkey>` (e.g. `\cite Schenke:2012wb`)
  inside an actual Doxygen doc comment (`/** ... */` or `///`), backed by
  `doc/ipglasma.bib`. Add a new entry to that file (INSPIRE-HEP's BibTeX
  export, e.g. `curl -H "Accept: application/x-bibtex"
  "https://inspirehep.net/api/arxiv/<id>"`) before citing a paper that
  isn't in it yet.
- `\cite` only has an effect inside a real Doxygen comment -- a citation
  in a plain `//` implementation comment (not extracted into the
  generated docs) should stay as plain text, e.g. `arXiv:1508.06294`.
- Top-level narrative docs (`README.md`, this file) intentionally keep
  plain Markdown links instead of `\cite`/`\iref`, since GitHub renders
  them directly and doesn't understand Doxygen's special commands
  (matching SMASH's own README).

### Scope

- Document every function, including trivial getters/setters -- don't
  default to grouping repetitive accessors under one shared comment.
- Document every class (its role in the simulation pipeline: who
  constructs it, what populates its state, what stage consumes it) and
  every member variable.
- Document every enum (what it selects/represents) and every enumerator
  individually, the same way as a struct and its member variables --
  don't leave an enum with only a plain `//` comment or no per-value
  docs.

### Example

```cpp
/**
 * Returns the projectile's (nucleus A) local \f$g^2\mu_A^2\f$ value.
 *
 * Used to sample the classical color source and to set the local
 * running coupling (RunningCoupling.cpp).
 * \return The stored \f$g^2\mu_A^2\f$ value [lattice units].
 */
double getg2mu2A() const { return g2mu2A_; }
```
