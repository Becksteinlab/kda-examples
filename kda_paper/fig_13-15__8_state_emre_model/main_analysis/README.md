## EmrE Main Analysis

Code to analyze the EmrE [free-exchange transport model](https://doi.org/10.1085/jgp.201912437) via [Kinetic Diagram Analysis](https://github.com/Becksteinlab/kda) (KDA).

**Warning:** KDA is in flux and not API stable. The latest known working KDA commit hash: `b63d5f2`


## Generating Figures

The modules here provide the ability to generate many different figures based on the 8-state model of EmrE. To generate _all_ figures, run the following:

```bash
$ mkdir plots
$ mkdir plots/figures
$ mkdir plots/flux_graphs
$ python plot_emre_8_state.py
$ python generate_figures.py
```

All generated figures can be found in either `plots/figures` or `plots/flux_graphs`. The raw data for each figure (7A, 7B, etc.) can be found in `data`. 

Here is a list of the key figures from the paper:

- `emre_8_state_model.png`: kinetic diagram for the 8-state model of EmrE
- `all_cycles.*`: All valid cycles for the 8-state model of EmrE (Figure 13)
- `fig_7A.*`: Operational fluxes of EmrE for alternating access rate biasing (Figure 14a)
- `fig_7B.*`: Operational fluxes of EmrE for substrate off-rate biasing (Figure 15a)
- `fig_7A_flux_diagram_RAA_1E-08.*`, `fig_7A_flux_diagram_RAA_1E+00.*`, `fig_7A_flux_diagram_RAA_1E+08.*`: kinetic diagrams with transition-flux-weighted edges for alternating access rate biasing (Figures 14b-14d). The legends for each are included separately (i.e. `fig_7A_flux_diagram_RAA_1E+00_legend.*`).
- `fig_7B_flux_diagram_kAA_1E+02_Roff_1E-10.*`, `fig_7B_flux_diagram_kAA_1E+02_Roff_1E+00.*`, `fig_7B_flux_diagram_kAA_1E+02_Roff_1E+02.*`: kinetic diagrams with transition-flux-weighted edges for substrate off-rate biasing (Figures 15b-15d). The legends for each are included separately.

**NOTE:** `generate_figures.py` uses previously generated data (`main_analysis/data/all_figure_data.csv`) and sympy expressions (located in `main_analysis/fig_functions`) to speed up code execution. `all_figure_data.csv` can be re-generated using `generate_fig_data.py`. To generate new sympy expressions at run time, simply delete all files in the `main_analysis/fig_functions` directory and re-run `generate_figures.py`. 
