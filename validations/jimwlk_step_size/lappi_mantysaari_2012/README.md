# Curves of arXiv:1212.4825

Curves of T. Lappi and H. Mäntysaari, "On the running coupling in the JIMWLK
equation", Eur. Phys. J. C 73 (2013) 2307,
[arXiv:1212.4825](https://arxiv.org/abs/1212.4825), used by
`../lappi_mantysaari.py`.

They were extracted from the vector graphics of the figure files in the
arXiv source: the vertices of each plotted line, converted to data
coordinates with the axis tick marks and their labels. Vertices added where
the plotting program clipped a line at the frame were dropped. The
square-root JIMWLK curve with c = 0.2 appears in both Fig. 8 and Fig. 11,
and the two extractions agree.

Each file has the columns `curve,x,y`:

| File | Figure (source file) | x | y | Curves |
|---|---|---|---|---|
| `fig6_dipole.csv` | Fig. 6 (`randrvssqrt_smallq_avgs25112012.pdf`) | rΛ_QCD | N = 1 − ⟨D(r)⟩ | `sqrtas_JIMWLK_y<y>` and `noiseas_JIMWLK_y<y>`, y = 0, 4, ..., 20 |
| `fig8_lambda.csv` | Fig. 8 (`lambdaranrdsrt_avgs25112012_wbk_wbc.pdf`) | Q_s/Λ_QCD | λ = d ln Q_s²/dy | JIMWLK with the square-root (`sqrtas`) and noise (`noiseas`) prescriptions and BK with the Balitsky prescription, for c = 0.2 (low and high initial Q_s) and c = 1.5 |
| `fig11_lambda.csv` | Fig. 11 (`lambdasqrtccomp.pdf`) | Q_s/Λ_QCD | λ | JIMWLK and BK with the square-root prescription, c = 0.2 and 1.5 |

The setup of the curves with the lower initial Q_s, which the script uses:
1024² lattice, the MV model with g²μL = 31 and N_y = 100, Λ_QCD L = 6,
μ₀L = 15, N_f = 3, no infrared regulator; Q_s from N(r_s) = 1 − e^{−1/2},
Q_s² = 2/r_s². The coupling has the scale 4e^{−2γ_E}/r².
