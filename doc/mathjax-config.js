// Extends Doxygen's default MathJax config (which only recognizes \(...\) --
// the form \f$...\f$ is converted to by Doxygen itself) to also recognize
// GitHub-flavored $...$ inline math, so README.md's own $...$ formulas (kept
// that way so GitHub renders them natively; see CONTRIBUTING.md) also render
// correctly here instead of showing up as raw, unprocessed LaTeX text.
MathJax.Hub.Config({
    tex2jax: {
        inlineMath: [
            ["$", "$"],
            ["\\(", "\\)"],
        ],
    },
});
