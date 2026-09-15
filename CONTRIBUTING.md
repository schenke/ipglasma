# Contributing

## Code documentation

IP-Glasma is starting to adopt [Doxygen](https://www.doxygen.nl/)-style
documentation, following the conventions used by
[SMASH](https://github.com/smash-transport/smash) (see its
[CONTRIBUTING.md](https://github.com/smash-transport/smash/blob/main/CONTRIBUTING.md)).
No `Doxyfile`/build target exists yet in this repo, but writing comments this
way now means they need no rework once one is added, and in the meantime
they still read as ordinary comments and IDE tooltips.

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

- Don't use `\iref`/`\cite` yet -- both require a bibliography file
  (Inspire-backed or not) that this repo doesn't have set up. Cite papers
  in plain text for now (matching how the existing prose comments already
  do it, e.g. `arXiv:1508.06294`), and switch the existing plain-text
  citations to `\iref`/`\cite` once a Doxyfile and `.bib` file exist.

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
