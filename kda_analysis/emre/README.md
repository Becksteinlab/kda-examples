## EmrE Analysis

Code to analyze the EmrE [free-exchange transport model](https://doi.org/10.1085/jgp.201912437) via [Kinetic Diagram Analysis](https://github.com/Becksteinlab/kda) (KDA).

## Usage

**Warning:** KDA is in flux and not API stable. The latest known working KDA commit hash: `16fcb78`

With KDA installed, the analysis code can be executed the following way:

```bash
# create directory structure for storing plots
mkdir plots
mkdir plots/figures
mkdir plots/flux_graphs
# generate all figures
python plot_emre_8_state.py
python generate_figures.py
```
