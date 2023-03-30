import os
import sys
import dill

dill.settings["recursive"] = True
import numpy as np
import pandas as pd
import networkx as nx
import matplotlib as mpl
from matplotlib.lines import Line2D
import matplotlib.pyplot as plt

from kda import graph_utils, calculations, plotting

import operational_flux
import plotting_emre
from fig_7a import data_7a, sub_dict_7a_values, plot_fig_7A
from fig_7b import data_7b, sub_dict_7b_values, plot_fig_7B
from fig_7d import data_7d, sub_dict_7d_values, plot_fig_7D
from fig_9 import data_9, sub_dict_9_values, plot_fig_9


def fmtscientific(x, sig_figs=0):
    """
    Formats a multiple of 10 in scientific notation (i.e. "$1 \times 10^3$").
    """
    # start with python scientific notation
    if sig_figs == 0:
        sci_str = f"{x:.0e}"
    elif sig_figs == 1:
        sci_str = f"{x:.1e}"
    # split the string into base and power
    base, pwr = sci_str.split("e")
    # cast as integer to remove the `+` signs and leading zeros
    pwr = str(int(pwr))
    # set in math text
    formatted_str = r"$" + base + r"\times 10^{" + pwr + r"}$"

    return formatted_str


def get_K(
    k31=1,
    k13=1,
    k57=1,
    k75=1,
    k42=1,
    k24=1,
    k68=1,
    k86=1,
    k34=1,
    k43=1,
    k56=1,
    k65=1,
    k12=1,
    k21=1,
    k78=1,
    k87=1,
    k71=1,
    k17=1,
    k53=1,
    k35=1,
    k64=1,
    k46=1,
    k82=1,
    k28=1,
):
    K = np.array(
        [
            [0, k12, k13, 0, 0, 0, k17, 0],
            [k21, 0, 0, k24, 0, 0, 0, k28],
            [k31, 0, 0, k34, k35, 0, 0, 0],
            [0, k42, k43, 0, 0, k46, 0, 0],
            [0, 0, k53, 0, 0, k56, k57, 0],
            [0, 0, 0, k64, k65, 0, 0, k68],
            [k71, 0, 0, 0, k75, 0, 0, k78],
            [0, k82, 0, 0, 0, k86, k87, 0],
        ]
    )
    return K


def get_K_and_G():
    # generate array using defaults (1)
    K = get_K()
    G = nx.MultiDiGraph()
    graph_utils.generate_edges(G, K)
    return K, G


def get_sub_dict(fig_key):
    if fig_key == "7A":
        return sub_dict_7a_values()
    elif fig_key == "7B":
        return sub_dict_7b_values()
    elif fig_key == "7D":
        return sub_dict_7d_values()
    elif fig_key == "9":
        return sub_dict_9_values()


def get_rate_names(fig_key):
    if fig_key == "7A":
        return ["k_AA_anti", "k_AA_sym"]
    elif fig_key == "7B":
        return ["k_EH_H", "k_EHD_H", "k_ED_D", "k_EHD_D", "k_conf"]
    elif fig_key == "7D":
        return ["k_EH_H", "k_EHD_H", "k_ED_D", "k_EHD_D", "k_AA_anti", "k_AA_sym"]
    elif fig_key == "9":
        return ["H_ext"]


def get_H_D_cos_dicts(G, K):
    # H: 1 -> 7, 2 -> 8
    # D: 4 -> 6, 2 -> 8
    H_edges = [(0, 6, 0), (1, 7, 0)]
    D_edges = [(3, 5, 0), (1, 7, 0)]

    H_cycles_valid = operational_flux.find_relevant_cycles(G, H_edges)[0]
    D_cycles_valid = operational_flux.find_relevant_cycles(G, D_edges)[0]
    H_cos_dict, D_cos_dict = operational_flux.get_cos_dicts(
        H_cycles_valid, D_cycles_valid
    )
    return H_cos_dict, D_cos_dict


def generate_sympy_func(G, cos_dict, sub_dict, filepath):
    sympy_func = operational_flux.get_op_cycle_flux_funcs(
        G=G, cos_dict=cos_dict, sub_dict=sub_dict, key="name"
    )
    dill.dump(sympy_func, open(filepath, "wb"))


def get_sympy_func(G, cos_dict, sub_dict, filepath):
    if not os.path.isfile(filepath):
        print(
            f"No SymPy function found at location {filepath} \n"
            f"Generating new sympy function..."
        )
        generate_sympy_func(
            G=G,
            cos_dict=cos_dict,
            sub_dict=sub_dict,
            filepath=filepath,
        )

    print(f"--> Loading SymPy function from location {filepath}")
    sympy_func = dill.load(open(filepath, "rb"))
    return sympy_func


def get_plot_data(df, fH, fD, fig_key):

    if fig_key == "7A":
        k_AA_anti_arr = df["k17"].values
        k_AA_sym_arr = df["k19"].values
        R_AA = k_AA_anti_arr / k_AA_sym_arr

        H_flux = np.zeros((R_AA.size), dtype=float)
        D_flux = np.zeros_like(H_flux)

        for i, (k_AA_anti, k_AA_sym) in enumerate(zip(k_AA_anti_arr, k_AA_sym_arr)):
            H_flux[i] = fH(k_AA_anti=k_AA_anti, k_AA_sym=k_AA_sym)
            D_flux[i] = fD(k_AA_anti=k_AA_anti, k_AA_sym=k_AA_sym)

        cols = ["R_AA", "Proton Flux (/s)", "Substrate Flux (/s)"]
        data = zip(R_AA, H_flux, D_flux)
        plot_df = pd.DataFrame(data=data, columns=cols)

    elif fig_key == "7B":

        k_EH_H_arr = df["k2"].values
        k_ED_D_arr = df["k10"].values
        k_EHD_H_arr = df["k6"].values
        k_EHD_D_arr = df["k14"].values
        k_AA_arr = df["k17"].values
        R_off = k_EH_H_arr / k_EHD_H_arr

        H_flux = np.zeros((R_off.size), dtype=float)
        D_flux = np.zeros_like(H_flux)
        for i, (k_EH_H, k_ED_D, k_EHD_H, k_EHD_D, k_conf) in enumerate(
            zip(k_EH_H_arr, k_ED_D_arr, k_EHD_H_arr, k_EHD_D_arr, k_AA_arr)
        ):
            H_flux[i] = fH(
                k_EH_H=k_EH_H,
                k_ED_D=k_ED_D,
                k_EHD_H=k_EHD_H,
                k_EHD_D=k_EHD_D,
                k_conf=k_conf,
            )
            D_flux[i] = fD(
                k_EH_H=k_EH_H,
                k_ED_D=k_ED_D,
                k_EHD_H=k_EHD_H,
                k_EHD_D=k_EHD_D,
                k_conf=k_conf,
            )

        cols = ["R_off", "k_AA", "Proton Flux (/s)", "Substrate Flux (/s)"]
        data = zip(R_off, k_AA_arr, H_flux, D_flux)
        plot_df = pd.DataFrame(data=data, columns=cols)

    elif fig_key == "7D":

        k_EH_H_arr = df["k2"].values
        k_ED_D_arr = df["k10"].values
        k_EHD_H_arr = df["k6"].values
        k_EHD_D_arr = df["k14"].values
        k_AA_anti_arr = df["k17"].values
        k_AA_sym_arr = df["k19"].values
        R_off = k_EH_H_arr / k_EHD_H_arr
        R_AA = k_AA_anti_arr / k_AA_sym_arr

        H_flux = np.zeros((R_off.size), dtype=float)
        D_flux = np.zeros_like(H_flux)
        for i, (k_EH_H, k_ED_D, k_EHD_H, k_EHD_D, k_AA_anti, k_AA_sym) in enumerate(
            zip(
                k_EH_H_arr,
                k_ED_D_arr,
                k_EHD_H_arr,
                k_EHD_D_arr,
                k_AA_anti_arr,
                k_AA_sym_arr,
            )
        ):
            H_flux[i] = fH(
                k_EH_H=k_EH_H,
                k_ED_D=k_ED_D,
                k_EHD_H=k_EHD_H,
                k_EHD_D=k_EHD_D,
                k_AA_anti=k_AA_anti,
                k_AA_sym=k_AA_sym,
            )
            D_flux[i] = fD(
                k_EH_H=k_EH_H,
                k_ED_D=k_ED_D,
                k_EHD_H=k_EHD_H,
                k_EHD_D=k_EHD_D,
                k_AA_anti=k_AA_anti,
                k_AA_sym=k_AA_sym,
            )

        cols = [
            "R_off",
            "k_AA_anti",
            "k_AA_sym",
            "R_AA",
            "Proton Flux (/s)",
            "Substrate Flux (/s)",
        ]
        data = zip(R_off, k_AA_anti_arr, k_AA_sym_arr, R_AA, H_flux, D_flux)
        plot_df = pd.DataFrame(data=data, columns=cols)

    elif fig_key == "9":

        cH_ext_arr = np.logspace(-5.4, -9.4, 209)
        pH_ext = -1 * np.log10(cH_ext_arr)

        H_flux = np.zeros(cH_ext_arr.size)
        D_flux = np.zeros_like(H_flux)
        for i, cH_ext in enumerate(cH_ext_arr):
            H_flux[i] = fH(H_ext=cH_ext)
            D_flux[i] = fD(H_ext=cH_ext)

        cols = ["External pH", "Proton Flux (/s)", "Substrate Flux (/s)"]
        data = zip(pH_ext, H_flux, D_flux)
        plot_df = pd.DataFrame(data=data, columns=cols)

    return plot_df


def plot_data(df, fig_key):

    colors = ["#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e"]

    if fig_key == "7A":
        fig, axs = plot_fig_7A(df=df)
    elif fig_key == "7B":
        fig, axs = plot_fig_7B(df=df, colors=colors)
    elif fig_key == "7D":
        fig, axs = plot_fig_7D(df=df, colors=colors)
    elif fig_key == "9":
        fig, axs = plot_fig_9(df=df)

    plot_stoichiometric_ratio(df=df, fig_key=fig_key, fig=fig, axs=axs)

    x_data_key = df.columns[0]
    if fig_key == "9":
        x_label = x_data_key
    else:
        x_label = r"$R_\mathrm{%s}$" % (x_data_key[2:])

    # mute the top subplot's xticklabels
    plt.setp(axs[-1].get_xticklabels(), visible=False)

    if fig_key == "7A":
        ax_leg = axs[1]
        ax_leg.axis("off")
        handles = []
        labels = []
        for _ax in fig.axes:
            _handles, _labels = _ax.get_legend_handles_labels()
            handles.extend(_handles)
            labels.extend(_labels)
        ax_leg.legend(handles, labels, loc="center", bbox_to_anchor=(0.5, 0.5))
        fig_leg = ax_leg.get_figure()
        for _ext in ("png", "pdf", "svg"):
            fig_leg.savefig(f"plots/figures/fig_7A_legend.{_ext}", dpi=300)

    elif fig_key == "9":
        fig.legend(loc="upper right", bbox_to_anchor=(0.97, 0.53))

    if fig_key in ["7A", "9"]:
        axs[0].set_ylabel(r"Flux (s$^{-1}$)")

    for ax in axs:
        ax.grid(True)
    return fig


def plot_stoichiometric_ratio(df, fig_key, fig, axs):
    colors = ["#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e"]

    if fig_key in ["7B", "7D"]:
        ax = fig.add_subplot(313, sharex=axs[0])
    else:
        ax = fig.add_subplot(212, sharex=axs[0])

    # ax.set_title(f"Fig. {fig_key} Stoichiometric Ratio", fontsize=11)

    column_keys = df.columns
    H_flux_key = column_keys[-2]
    D_flux_key = column_keys[-1]

    H_str = r"$J_{\mathrm{H}^{+}}$"
    D_str = r"$J_\mathrm{D}$"
    label = D_str + "/" + H_str

    if fig_key == "7A":
        x_vals = df["R_AA"].values
        y_vals = df[D_flux_key].values / df[H_flux_key].values
        ax.plot(x_vals, y_vals, color="dimgrey", label=label)
        ax.axvline(
            x=1,
            ymin=0,
            ymax=1,
            ls="--",
            color="black",
            # label=r"$R_\mathrm{AA}$ = 1",
        )
        legend_title = None

    elif fig_key == "7B":
        x_vals = df["R_off"].values
        y_vals = df[D_flux_key].values / df[H_flux_key].values
        k_AA_data = df["k_AA"].values
        unique_k_AA = np.unique(k_AA_data)

        for i, val in enumerate(unique_k_AA):
            mask = k_AA_data == val
            label = f"{val:.0e}"
            ax.plot(
                x_vals[mask],
                y_vals[mask],
                ls="-",
                lw=0.6,
                color=colors[i],
                label=label,
            )
        ax.axvline(
            x=1,
            ymin=0,
            ymax=1,
            ls="--",
            lw=0.8,
            color="black",
            # label=r"$R_\mathrm{off}$ = 1",
        )
        legend_title = r"$k_\mathrm{AA}$"

    elif fig_key == "7D":
        x_vals = df["R_off"].values
        y_vals = df[D_flux_key].values / df[H_flux_key].values
        R_AA_data = df["R_AA"].values
        unique_R_AA = np.unique(R_AA_data)

        for i, val in enumerate(unique_R_AA):
            mask = R_AA_data == val
            label = f"{val:.1f}"
            ax.plot(
                x_vals[mask],
                y_vals[mask],
                ls="-",
                lw=0.6,
                color=colors[i],
                label=label,
            )
        ax.axvline(
            x=1,
            ymin=0,
            ymax=1,
            ls="--",
            lw=0.8,
            color="black",
            label=r"$R_\mathrm{off}$ = 1",
        )
        legend_title = r"$R_\mathrm{AA}$"
    elif fig_key == "9":
        pH_key = df.columns[0]
        x_vals = df[pH_key].values
        y_vals = df[D_flux_key].values / df[H_flux_key].values
        ax.semilogx(x_vals, y_vals, color="dimgrey", label=label)
        ax.axvline(x=7.4, ymin=0, ymax=1, ls="--", color="black", label=r"pH = 7.4")
        x_ticks = np.linspace(5.4, 9.4, 5)
        ax.set_xticks(x_ticks, minor=False)
        ax.set_xticks([], minor=True)
        ax.set_xticklabels(x_ticks)
        legend_title = None

    x_data_key = df.columns[0]
    if fig_key == "9":
        x_label = x_data_key
    else:
        x_label = r"$R_\mathrm{%s}$" % (x_data_key[2:])
    ax.set_xlabel(x_label)
    ax.set_ylabel("Stoichiometric Ratio")
    ax.grid(True)

    if fig_key != "9":
        ax.set_xscale("log")
        ax.set_xticks([])
        ax.set_xticklabels([])
        max_power = np.max(np.log10(x_vals))
        if fig_key == "7A":
            xticks = np.logspace(-max_power, max_power, 7)
        else:
            xticks = np.logspace(-max_power, max_power, 9)
        ax.set_xticks(xticks)
        ax.set_xticklabels(xticks)
        logfmt = mpl.ticker.LogFormatterSciNotation(base=10.0, labelOnlyBase=True)
        ax.xaxis.set_major_formatter(logfmt)


def get_unique_idx_pairs(G):
    indices = []
    for i, j, k in G.edges:
        if (j, i, k) not in indices:
            indices.append((i, j, k))
    return indices


def get_trans_flux_range(G, K, state_probs):

    trans_fluxes = []
    net_trans_fluxes = []
    for i, j, k in G.edges:
        J_for = state_probs[i] * K[i, j]
        J_rev = state_probs[j] * K[j, i]
        net_trans_flux = J_for - J_rev

        trans_fluxes.append(J_for)
        trans_fluxes.append(J_rev)
        if np.abs(net_trans_flux) > 1e-15:
            if net_trans_flux < 0:
                net_trans_fluxes.append(-net_trans_flux)
            else:
                net_trans_fluxes.append(net_trans_flux)

    min_tflux = np.min(trans_fluxes)
    max_tflux = np.max(trans_fluxes)
    min_ntflux = np.min(net_trans_fluxes)
    max_ntflux = np.max(net_trans_fluxes)

    trans_flux_range = (min_tflux, max_tflux)
    net_trans_flux_range = (min_ntflux, max_ntflux)
    return trans_flux_range, net_trans_flux_range


def get_minmax_data(x_array, xmin, xmax, a=0.6, b=5, logscale=False):
    """
    Smallest edge width (a): 0.6
    Largest edge width (b): 5

    f(xmin) = a
    f(xmax) = b

           (b - a)(x - xmin)
    f(x) = -----------------  + a
              xmax - xmin
    """
    if logscale:
        x_array = np.log10(x_array)
        xmin = np.log10(xmin)
        xmax = np.log10(xmax)
    minmax_data = np.zeros_like(x_array)
    for i, x in enumerate(x_array):
        numerator = (b - a) * (x - xmin)
        denominator = xmax - xmin
        val = (numerator / denominator) + a
        minmax_data[i] = val
    return minmax_data


def generate_flux_graph(
    G,
    K,
    probs,
    tflux_range,
    ntflux_range,
    title,
    fig_key,
    node_pos,
):
    trans_fluxes = []
    trans_flux_edges = []
    for (i, j, k) in G.edges:
        J = probs[i] * K[i, j]
        trans_fluxes.append(J)
        trans_flux_edges.append((i, j))

    unique_idx = get_unique_idx_pairs(G=G)

    net_trans_fluxes = []
    net_trans_flux_edges = []
    ntflux_label_dict = {}
    for (i, j, k) in unique_idx:
        J_for = probs[i] * K[i, j]
        J_rev = probs[j] * K[j, i]
        net_trans_flux = J_for - J_rev

        if np.abs(net_trans_flux) > 0.02:
            val = np.abs(net_trans_flux)
            label = f"{val:.2f}"
            ntflux_label_dict[(i, j)] = label

        if np.abs(net_trans_flux) > 1e-15:
            if net_trans_flux < 0:
                net_trans_fluxes.append(-net_trans_flux)
                net_trans_flux_edges.append((j, i))
            else:
                net_trans_fluxes.append(net_trans_flux)
                net_trans_flux_edges.append((i, j))

    scale_factor = 1.25
    diagram_side = 1.625 * scale_factor # inches
    legend_height = 0.875 * scale_factor # inches
    legend_width = 1.25 * scale_factor # inches

    fig = plt.figure(figsize=(diagram_side, diagram_side))
    ax = fig.add_subplot(111)
    ax.axis("off")

    fig_leg = plt.figure(figsize=(legend_width, legend_height))
    ax_leg = fig_leg.add_subplot(111)
    ax_leg.axis("off")

    node_list = list(G.nodes)

    arrowsize = 8
    arrowstyle = "->"
    arrow_color = "0.8"
    ntflux_color = "red"

    node_size = get_minmax_data(x_array=probs, xmin=0, xmax=1, a=60, b=500)
    tflux_weights = get_minmax_data(
        x_array=trans_fluxes,
        xmin=tflux_range[0],
        xmax=tflux_range[1],
        a=1.5,
        b=3,
        logscale=True,
    )
    ntflux_weights = get_minmax_data(
        x_array=net_trans_fluxes,
        xmin=ntflux_range[0],
        xmax=ntflux_range[1],
        a=0.2,
        b=1.5,
        logscale=False,
    )

    node_labels = plotting._get_node_labels(node_list=node_list)
    bbox = dict(boxstyle="square, pad=0.2", fc="white", ec="black", lw=0.5, alpha=0.7)

    nx.draw_networkx_nodes(
        G, node_pos, node_size=node_size, nodelist=node_list, node_color=arrow_color, ax=ax
    )

    nx.draw_networkx_edges(
        G,
        node_pos,
        edgelist=trans_flux_edges,
        node_size=node_size * 2,
        width=tflux_weights,
        edge_color=arrow_color,
        arrowsize=arrowsize,
        arrowstyle=arrowstyle,
        connectionstyle="arc3, rad = 0.11",
        ax=ax,
    )

    nx.draw_networkx_edges(
        G,
        node_pos,
        edgelist=net_trans_flux_edges,
        node_size=node_size * 1.5,
        width=ntflux_weights,
        edge_color=ntflux_color,
        arrowsize=arrowsize,
        arrowstyle=arrowstyle,
        ax=ax,
    )

    nx.draw_networkx_labels(
        G,
        node_pos,
        node_labels,
        font_size=7,
        horizontalalignment="center",
        verticalalignment="center",
        ax=ax,
    )

    nx.draw_networkx_edge_labels(
        G,
        pos=node_pos,
        edge_labels=ntflux_label_dict,
        font_size=5.5,
        font_color="black",
        bbox=bbox,
        ax=ax,
    )

    J_min = r"$J_\mathrm{min} =$"
    J_max = r"$J_\mathrm{max} =$"
    tflux_min_label = J_min + fmtscientific(np.min(trans_fluxes), sig_figs=1)
    tflux_max_label = J_max + fmtscientific(np.max(trans_fluxes), sig_figs=1)
    ntflux_min_label = r"$\Delta$ " + J_min + fmtscientific(np.min(net_trans_fluxes), sig_figs=1)
    ntflux_max_label = r"$\Delta$ " + J_max + fmtscientific(np.max(net_trans_fluxes), sig_figs=1)

    legend_elements = [
        Line2D(
            [0],
            [0],
            color=arrow_color,
            lw=tflux_weights.min(),
            label=tflux_min_label,
        ),
        Line2D(
            [0],
            [0],
            color=arrow_color,
            lw=tflux_weights.max(),
            label=tflux_max_label,
        ),
        Line2D(
            [0],
            [0],
            color=ntflux_color,
            lw=ntflux_weights.min(),
            label=ntflux_min_label,
        ),
        Line2D(
            [0],
            [0],
            color=ntflux_color,
            lw=ntflux_weights.max(),
            label=ntflux_max_label,
        ),
    ]

    ax_leg.legend(
        loc="center",
        bbox_to_anchor=(0.49, 0.5),
        handles=legend_elements,
        title=title,
    )
    fig.subplots_adjust(bottom=0, top=1, left=0, right=1)
    return fig, fig_leg


def plot_flux_graphs(df, fig_key):

    if fig_key == "7A":
        n_datasets = df.shape[0]
        subset_idx = np.linspace(0, n_datasets - 1, 5, dtype=np.int32)
        df_subset = df.iloc[subset_idx]
    elif fig_key == "7B":
        # for this figure there are several values of k_AA to choose from,
        # but we don't really need to compare between AA rates, so let's
        # just select one of the middle values
        selected_k_AA = 100
        new_df = df[df["k17"] == selected_k_AA]
        R_off = new_df["k2"] / new_df["k6"]
        R_off = R_off.reset_index(drop=True)
        # values of interest
        interest_vals = [1e-10, 1, 100]
        # find the indices of the dataframe where R_off values
        # are equal to the values of interest
        subset_idx = [R_off.index[R_off == val].tolist()[0] for val in interest_vals]
        df_subset = new_df.iloc[subset_idx]
    else:
        raise NotImplementedError(f"Not implemented for Fig. {fig_key}")

    node_pos = plotting_emre.get_node_positions()

    H_in = 10 ** (-6.5)
    H_out = 10 ** (-7.5)
    c_D = 25e-9

    data_dict = {}
    data_dict["graphs"] = []
    data_dict["rate_matrices"] = []
    data_dict["fig_titles"] = []
    data_dict["filenames"] = []
    data_dict["probabilities"] = []
    data_dict["trans_flux_ranges"] = []
    data_dict["net_trans_flux_ranges"] = []
    data_dict["max_flux_case"] = []
    for idx, row in df_subset.iterrows():
        K = get_K(
            row["k1"] * H_out,
            row["k2"],
            row["k3"] * H_in,
            row["k4"],
            row["k5"] * H_out,
            row["k6"],
            row["k7"] * H_in,
            row["k8"],
            row["k9"] * c_D,
            row["k10"],
            row["k11"] * c_D,
            row["k12"],
            row["k13"] * c_D,
            row["k14"],
            row["k15"] * c_D,
            row["k16"],
            row["k17"],
            row["k18"],
            row["k19"],
            row["k20"],
            row["k21"],
            row["k22"],
            row["k23"],
            row["k24"],
        )
        G = nx.MultiDiGraph()
        graph_utils.generate_edges(G, K)
        state_probs = calculations.calc_state_probs(G, key="val")
        trans_flux_range, net_trans_flux_range = get_trans_flux_range(
            G=G, K=K, state_probs=state_probs
        )
        data_dict["graphs"].append(G)
        data_dict["rate_matrices"].append(K)
        data_dict["probabilities"].append(state_probs)
        data_dict["trans_flux_ranges"].append(trans_flux_range)
        data_dict["net_trans_flux_ranges"].append(net_trans_flux_range)

        is_max_flux = False

        if fig_key == "7A":
            R_AA = row["k17"] / row["k19"]
            if R_AA == 1:
                is_max_flux = True
            title_str = r"$R_\mathrm{AA} = $" + fmtscientific(R_AA)
            filename = f"fig_{fig_key}_flux_diagram_RAA_{R_AA:.0E}"
            filename_leg = f"fig_{fig_key}_flux_diagram_RAA_{R_AA:.0E}_legend"
        elif fig_key == "7B":
            R_off = row["k2"] / row["k6"]
            if R_off == 1:
                is_max_flux = True
            k_AA = row["k17"]
            title_str = r"$R_\mathrm{off} = $" + fmtscientific(R_off)
            filename = f"fig_{fig_key}_flux_diagram_kAA_{k_AA:.0E}_Roff_{R_off:.0E}"
            filename_leg = f"fig_{fig_key}_flux_diagram_kAA_{k_AA:.0E}_Roff_{R_off:.0E}_legend"

        data_dict["max_flux_case"].append(is_max_flux)
        data_dict["filenames"].append(filename)
        data_dict["filenames"].append(filename_leg)
        data_dict["fig_titles"].append(title_str)

    sorted_tflux = np.sort(np.array(data_dict["trans_flux_ranges"]).flatten())
    sorted_ntflux = np.sort(np.array(data_dict["net_trans_flux_ranges"]).flatten())

    data_dict["figures"] = []
    for (G, K, probs, title, is_max_flux_case) in zip(
        data_dict["graphs"],
        data_dict["rate_matrices"],
        data_dict["probabilities"],
        data_dict["fig_titles"],
        data_dict["max_flux_case"],
    ):
        if is_max_flux_case:
            tflux_range = (sorted_tflux[0], sorted_tflux[-1])
            ntflux_range = (sorted_ntflux[0], sorted_ntflux[-1])
        else:
            max_flux_mask = ~np.array(data_dict["max_flux_case"])
            tmp_tflux_data = np.array(data_dict["trans_flux_ranges"])[max_flux_mask]
            tmp_ntflux_data = np.array(data_dict["net_trans_flux_ranges"])[
                max_flux_mask
            ]
            tmp_sorted_tflux = np.sort(tmp_tflux_data.flatten())
            tmp_sorted_ntflux = np.sort(tmp_ntflux_data.flatten())
            tflux_range = (tmp_sorted_tflux[0], tmp_sorted_tflux[-1])
            ntflux_range = (tmp_sorted_ntflux[0], tmp_sorted_ntflux[-1])

        fig, fig_leg = generate_flux_graph(
            G=G,
            K=K,
            probs=probs,
            tflux_range=tflux_range,
            ntflux_range=ntflux_range,
            title=title,
            fig_key=fig_key,
            node_pos=node_pos,
        )
        data_dict["figures"].append(fig)
        data_dict["figures"].append(fig_leg)
        plt.close()

    return data_dict
