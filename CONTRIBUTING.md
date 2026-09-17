# Contributing

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
  docs (a gap this pass initially left on `IntegrandId`/`NucleusRole`
  in Glauber.h, added before this rule was; caught by the user).

### Example

```cpp
/**
 * Local \f$g^2\mu_A^2\f$ color-charge-density-squared value for the
 * projectile (nucleus A) at this site, sampled during initialization.
 *
 * Feeds the classical color source and the local running-coupling/
 * saturation-scale determination in Evolution.cpp.
 * \return The stored \f$g^2\mu_A^2\f$ value [lattice units].
 */
double getg2mu2A() const { return g2mu2A_; }
```
