import os
import numpy as np
import networkx as nx

from kda import (
    diagrams,
    calculations,
    graph_utils,
    plotting as plotting_kda,
)

import fig_7a
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

    ##################
    ##### Plot G #####
    ##################
    plotting_kda.draw_diagrams(G, node_pos, path=save_path, label="emre_8_state_model")

    #######################
    ##### Plot Cycles #####
    #######################
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

    #################################
    ##### Plot Net Cycle Fluxes #####
    #################################

    dir_path = "fig_functions"
    test_path = os.path.join(dir_path, "fig_7A_cycle_1_net_cycle_flux.txt")
    if not os.path.isfile(test_path):
        sub_dict = fig_7a.sub_dict_7a_values()
        op_flux.generate_net_cycle_flux_sympy_funcs(G, sub_dict=sub_dict, dir_path=dir_path)

    df = fig_7a.get_7a_dataframe()
    k_AA_anti_arr = df["k17"].values
    k_AA_sym_arr = df["k19"].values
    R_AA_arr = k_AA_anti_arr / k_AA_sym_arr

    cycle_data_dict = op_flux.get_net_cycle_flux_cycles_and_funcs(
        G, rate_names=["k_AA_anti", "k_AA_sym"], dir_path=dir_path
    )

    flux_arrays = []
    for cycle_key, data_dict in cycle_data_dict.items():
        cycle = data_dict["cycle"]
        func = data_dict["func"]
        flux_arr = np.asarray(list(map(func, k_AA_anti_arr, k_AA_sym_arr)), dtype=np.float64)
        flux_arrays.append(flux_arr)

    flux_total = np.asarray(flux_arrays).sum(axis=0)

    # plot the individual net cycle fluxes
    cycle_flux_fig = plotting_emre.plot_net_cycle_fluxes(
        R=R_AA_arr,
        fluxes=flux_arrays,
        flux_tot=flux_total,
        cycles=all_cycles,
        norm=True,
    )
    for ext in ("png", "svg", "pdf"):
        cycle_flux_fig.savefig(os.path.join(save_path, f"fig_7A_net_cycle_fluxes.{ext}"), dpi=300)
