import os
import networkx as nx

from kda import (
    graph_utils,
    plotting as plotting_kda,
)

import fig_utils
import operational_flux as op_flux
import plotting_emre

if __name__ == "__main__":

    K = fig_utils.get_K()
    G = nx.MultiDiGraph()
    graph_utils.generate_edges(G, K)

    all_cycles = graph_utils.find_all_unique_cycles(G)
    H_edges = [(0, 6, 0), (1, 7, 0)]
    D_edges = [(3, 5, 0), (1, 7, 0)]
    H_cycles_valid, H_cycles_invalid = op_flux.find_relevant_cycles(G, H_edges)
    D_cycles_valid, D_cycles_invalid = op_flux.find_relevant_cycles(G, D_edges)

    node_pos = plotting_emre.get_node_positions()

    # make sure to run it from ~/EmrE directory
    save_path = os.path.join(os.getcwd(), "plots/figures/")
    print(f"Saving plots at location: {save_path}")
    plotting_kda.draw_diagrams(G, node_pos, path=save_path, label="emre_8_state_model")

    plotting_kda.draw_cycles(
        G,
        all_cycles,
        pos=node_pos,
        panel=True,
        panel_scale=2,
        font_size=13,
        rows=4,
        cbt=True,
        curved_arrows=True,
        path=save_path,
        label="all_cycles",
    )

    plotting_kda.draw_cycles(
        G,
        H_cycles_valid,
        pos=node_pos,
        panel=True,
        panel_scale=2,
        cbt=True,
        path=save_path,
        label="H_cycles_valid",
    )
    plotting_kda.draw_cycles(
        G,
        D_cycles_valid,
        pos=node_pos,
        panel=True,
        panel_scale=2,
        cbt=True,
        path=save_path,
        label="D_cycles_valid",
    )
    plotting_kda.draw_cycles(
        G,
        H_cycles_invalid,
        pos=node_pos,
        panel=True,
        panel_scale=2,
        cbt=True,
        path=save_path,
        label="H_cycles_invalid",
    )
    plotting_kda.draw_cycles(
        G,
        D_cycles_invalid,
        pos=node_pos,
        panel=True,
        panel_scale=2,
        cbt=True,
        path=save_path,
        label="D_cycles_invalid",
    )
