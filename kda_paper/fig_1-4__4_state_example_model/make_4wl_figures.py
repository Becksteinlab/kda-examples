# generate 4wl diagrams

import os
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt

from kda import plotting, graph_utils, diagrams


def _flatten_flux_diagrams(diagrams):
    # temporary work around to flatten output of
    # `kda.diagrams.generate_all_flux_diagrams()`
    flux_diagrams = []
    for diag_list in diagrams:
        if diag_list:
            for diag in diag_list:
                flux_diagrams.append(diag)
    return flux_diagrams


def save_figure(fig, filename):
    # get the current working directory for saving figures
    cwd = os.getcwd()
    for _ext in ("png", "svg", "pdf"):
        path = os.path.join(cwd, f"{filename}.{_ext}")
        fig.savefig(path, dpi=300)


def main():
	K = np.array([
	    [0, 1, 0, 1],
	    [1, 0, 1, 1],
	    [0, 1, 0, 1],
	    [1, 1, 1, 0],
	])

	G = nx.MultiDiGraph()
	graph_utils.generate_edges(G, K)
	pos = {0 : [1, 1],
	       1 : [-1, 1],
	       2 : [-1, -1],
	       3 : [1, -1]}

	pars = diagrams.generate_partial_diagrams(G)
	dirpars = diagrams.generate_directional_diagrams(G)
	flux_diagrams = diagrams.generate_all_flux_diagrams(G)
	flux_diagrams = _flatten_flux_diagrams(flux_diagrams)

	G_fig = plotting.draw_diagrams(G, pos=pos, font_size=12, curved_arrows=False)
	pars_fig = plotting.draw_diagrams(pars, pos=pos, panel=True, panel_scale=1.75, rows=2, font_size=12)
	dirs_fig = plotting.draw_diagrams(dirpars, pos=pos, panel=True, panel_scale=1.75, rows=4, font_size=12, cbt=True)
	flux_fig = plotting.draw_diagrams(
    	flux_diagrams, pos=pos, panel=True, panel_scale=1.75, 
    	rows=2, font_size=12, cbt=True, curved_arrows=True)

	figs = [G_fig, pars_fig, dirs_fig, flux_fig]
	fnames = ["4wl_diagram", "4wl_partial_diagram_panel", "4wl_directional_diagram_panel", "4wl_flux_diagram_panel"]
	for fig, fname in zip(figs, fnames):
	    save_figure(fig, fname)

if __name__ == "__main__":
	main()