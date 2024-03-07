# 6 State Model Net Cycle Flux Testing

import os
import numpy as np
import pandas as pd
import networkx as nx

import matplotlib.pyplot as plt
import sympy

from kda import plotting, graph_utils, calculations, diagrams, expressions, ode


def get_G1_rate_matrix():
    K = np.array(
        [
            [0, 1, 0, 0, 0, 1],
            [1, 0, 1, 0, 0, 0],
            [0, 1, 0, 1, 0, 0],
            [0, 0, 1, 0, 1, 0],
            [0, 0, 0, 1, 0, 1],
            [1, 0, 0, 0, 1, 0],
        ]
    )
    return K


def get_G2_rate_matrix():
    K = np.array(
        [
            [0, 1, 0, 1, 0, 1],
            [1, 0, 1, 0, 0, 0],
            [0, 1, 0, 1, 0, 0],
            [1, 0, 1, 0, 1, 0],
            [0, 0, 0, 1, 0, 1],
            [1, 0, 0, 0, 1, 0],
        ]
    )
    return K


def get_node_positions():
    pos = {}
    offset = (np.pi / 2)
    for i in range(6):
        theta = offset - i * (np.pi / 3)
        pos[i] = [np.cos(theta), np.sin(theta)]
    return pos


def generate_input_diagram(K):
    G = nx.MultiDiGraph()
    graph_utils.generate_edges(G, K)
    return G

def get_G1_substitution_dict():
    rate_names = ["A_on, A_off, B_on, B_off, k_AA_H, k_AA_Na, A_in, A_out, B_in, B_out"]
    A_on, A_off, B_on, B_off, k_AA_H, k_AA_Na, A_in, A_out, B_in, B_out = sympy.symbols('A_on A_off B_on B_off k_AA_H k_AA_Na A_in A_out B_in B_out')
    k12, k21, k23, k32, k34, k43, k45, k54, k56, k65, k61, k16 = sympy.symbols("k12 k21 k23 k32 k34 k43 k45 k54 k56 k65 k61 k16")
    sub_dict = {
        k12 : A_on*A_out,
        k21 : A_off,
        k23 : k_AA_H,
        k32 : k_AA_H,
        k34 : A_off,
        k43 : A_on*A_in,
        k45 : B_on*B_in,
        k54 : B_off,
        k56 : k_AA_Na,
        k65 : k_AA_Na,
        k61 : B_off,
        k16 : B_on*B_out,
    }
    return rate_names, sub_dict


def get_G2_substitution_dict():
    rate_names = ["A_on, A_off, B_on, B_off, k_AA_H, k_AA_Na, k_leak, A_in, A_out, B_in, B_out"]
    A_on, A_off, B_on, B_off, k_AA_H, k_AA_Na, k_leak,  A_in, A_out, B_in, B_out = sympy.symbols('A_on A_off B_on B_off k_AA_H k_AA_Na k_leak A_in A_out B_in B_out')
    k12, k21, k23, k32, k34, k43, k45, k54, k56, k65, k61, k16, k14, k41 = sympy.symbols("k12 k21 k23 k32 k34 k43 k45 k54 k56 k65 k61 k16 k14, k41")
    sub_dict = {
        k12 : A_on*A_out,
        k21 : A_off,
        k23 : k_AA_H,
        k32 : k_AA_H,
        k34 : A_off,
        k43 : A_on*A_in,
        k45 : B_on*B_in,
        k54 : B_off,
        k56 : k_AA_Na,
        k65 : k_AA_Na,
        k61 : B_off,
        k16 : B_on*B_out,
        k14 : k_leak,
        k41 : k_leak,
    }
    return rate_names, sub_dict


def get_G1_net_cycle_flux_lambda_function(G1):
    G1_order = [0, 1]
    G1_cycle = graph_utils.find_all_unique_cycles(G1)[0]
    G1_cycle_flux_str = calculations.calc_net_cycle_flux(G1, G1_cycle, order=G1_order, key='name', output_strings=True)
    rate_names1, sub_dict1 = get_G1_substitution_dict()
    J_sympy = G1_cycle_flux_str.subs(sub_dict1).simplify()
    J_lambda = expressions.construct_lambda_funcs(J_sympy, rate_names1)
    return J_lambda


def get_G2_net_cycle_flux_lambda_function(G2, G2_cycles):
    G2_orders = [[5, 0], [0, 1], [0, 1]]
    rate_names2, sub_dict2 = get_G2_substitution_dict()

    G2_J_sympy = []
    for cycle, order in zip(G2_cycles_reordered, G2_orders):
        J_str = calculations.calc_net_cycle_flux(G2, cycle, order=order, key='name', output_strings=True)
        J_sympy = J_str.subs(sub_dict2).simplify()
        G2_J_sympy.append(J_sympy)

    G2_J_lambda = expressions.construct_lambda_funcs(G2_J_sympy, rate_names2)
    return G2_J_lambda


def draw_cycles_custom(
    G,
    cycles,
    pos=None,
    panel_scale=2,
    rows=None,
    cols=None,
    font_size=12,
    cbt=False,
    curved_arrows=False):
    """
    Custom cycle plotting function based on the `kda.plotting.draw_cycles` function.
    """

    if curved_arrows:
        connection_style = "arc3, rad = 0.11"
    else:
        connection_style = "arc3"

    nrows, ncols, excess_plots = plotting._get_panel_dimensions(n_diagrams=len(cycles), rows=rows, cols=cols)

    fig, ax = plt.subplots(nrows=nrows, ncols=ncols, tight_layout=True)
    fig.set_figheight(nrows * panel_scale)
    fig.set_figwidth(1.2 * ncols * panel_scale)

    xlims, ylims = plotting._get_axis_limits(pos, scale_factor=1.4)

    for i, cycle in enumerate(cycles):
        node_labels = plotting._get_node_labels(node_list=cycle)
        node_colors = plotting._get_node_colors(cbt=cbt, obj=cycle)
        cycle_edges = plotting._construct_cycle_edges(cycle)
        edge_list = plotting._append_reverse_edges(cycle_edges)

        ix = np.unravel_index(i, ax.shape)
        plt.sca(ax[ix])
        ax[ix].set_xlim(xlims)
        ax[ix].set_ylim(ylims)

        plotting._plot_single_diagram(
            diagram=G,
            pos=pos,
            edge_list=edge_list,
            node_list=cycle,
            node_labels=node_labels,
            node_colors=node_colors,
            node_size=150 * panel_scale,
            font_size=font_size,
            arrow_width=1.5,
            cbt=False,
            connection_style=connection_style,
            ax=ax[ix],
        )

    for j in range(excess_plots):
        ax.flat[-j - 1].set_visible(False)
    return fig


def plot_stacked_flux_fig(leak_arr, G1_flux_arr, G2_productive_cycle_flux_arr, A_flux_arr, B_flux_arr):
    """
    Custom plotting function for plotting the net cycle flux,
    operational fluxes, and stoichiometry for both G1 and G2.
    """

    # colorblind palette
    prod_cycle_color = "#1b9e77"
    H_color = "#d95f02"
    Na_color = "#7570b3"

    fig, axs = plt.subplots(3, figsize=(5, 5), sharex=True, tight_layout=True)
    # calculate the transport stoichiometry
    stoich_arr = B_flux_arr/A_flux_arr
    # set the min/max values for the operational flux figures
    flux_figs_miny = -50
    flux_figs_maxy = 50

    print_vals = True
    if print_vals:
        print(f"G net cycle flux: {G1_flux_arr[0]:.3f}")
        print(f"G_leak net cycle flux (productive cycle) final: {G2_productive_cycle_flux_arr[-1]:.3f}")
        print(f"H+ operational flux final: {A_flux_arr[-1]:.3f}")
        print(f"Na+ operational flux final: {B_flux_arr[-1]:.3f}")
        print(f"Stoichiometry final: {stoich_arr[-1]:.3f}")

    axs[0].semilogx(leak_arr, G1_flux_arr, '--', lw=2, color="black", label="G")
    axs[0].semilogx(leak_arr, G2_productive_cycle_flux_arr, '-', lw=2, color=prod_cycle_color, label=r"G$_{\mathrm{leak}}$")
    axs[0].set_title("Net Cycle Flux: Coupled Cycle")
    axs[0].set_ylabel(r"Flux ($s^{-1}$)")
    axs[0].set_ylim(flux_figs_miny, flux_figs_maxy)
    axs[0].legend(loc="upper left", bbox_to_anchor=(1, 1))
    axs[0].grid(True)

    axs[1].semilogx(leak_arr, G1_flux_arr, '--', lw=2, color="black", label="G")
    axs[1].semilogx(leak_arr, A_flux_arr, '-', color=H_color, lw=2, label=r"$\mathrm{H}^+$")
    axs[1].semilogx(leak_arr, B_flux_arr, '-', color=Na_color, lw=2, label=r"$\mathrm{Na}^+$")
    axs[1].set_title("Operational Flux")
    axs[1].set_ylabel(r"Flux ($s^{-1}$)")
    axs[1].set_ylim(flux_figs_miny, flux_figs_maxy)
    axs[1].legend(loc="upper left", bbox_to_anchor=(1, 1))
    axs[1].grid(True)

    axs[2].axhline(y=1, ls='--', lw=2, color='black', label="G")
    axs[2].semilogx(leak_arr, stoich_arr, '-', lw=2, color="dimgrey", label=r"G$_{\mathrm{leak}}$")
    axs[2].set_title(r"Stoichiometry: Sodium Ion/Proton")
    axs[2].set_xlabel(r"$k_{\mathrm{leak}}$ ($s^{-1}$)")
    axs[2].set_ylim(-1.1, 1.1)
    axs[2].legend(loc="upper left", bbox_to_anchor=(1, 1))
    axs[2].grid(True)

    axs[0].set_xlim(leak_arr[0], leak_arr[-1])
    return fig


def create_antiporter_df(filename="antiporter_data.csv"):
    # constants
    NA = 6.022e23   # Avogadro's number
    D_H = 9.3e11    # Angstrom ^ 2 / s, diffusion coefficient for H+
    D_Na = 2.5e11   # Angstrom ^ 2 / s, diffusion coefficient for Na+
    R = 5           # Angstrom, approximate value using VMD
    A_conv = 1e-27  # Conversion from cubic angstroms (R*D) to liters

    # concentrations
    pH_int = 8.5
    pH_ext = 5.5
    # H_int < H_ext
    c_H_int = 10 ** (-pH_int)  # M
    c_H_ext = 10 ** (-pH_ext)  # M
    assert c_H_int < c_H_ext
    # Na_int < Na_ext
    c_Na_int = 10 * 1e-3    # M
    c_Na_ext = 150 * 1e-3   # M
    assert c_Na_int < c_Na_ext

    # conformational change rate (same both directions)
    k_AA_H = 69     # s^-1
    k_AA_Na = 350   # s^-1

    # binding/unbinding rates for H+

    flux = 4 * R * D_H                         # Weber disk solution
    k_on_prime_H = flux * A_conv * NA          # M^-1 s^-1
    pKa = 6.8
    k_off_H = k_on_prime_H * (10 ** (-pKa))     # s^-1

    # binding/unbinding rates for Na+
    Kd = 30 * 1e-3      # M, dissociation constant for Na+ binding

    flux = 4 * R * D_Na                         # Weber disk solution
    k_on_prime_Na = flux * A_conv * NA          # M^-1 s^-1
    k_off_Na = k_on_prime_Na * Kd               # s^-1

    data_dict = {
        "H_on": k_on_prime_H,
        "H_off": k_off_H,
        "Na_on": k_on_prime_Na,
        "Na_off": k_off_Na,
        "k_AA_H": k_AA_H,
        "k_AA_Na": k_AA_Na,
        "H_in": c_H_int,
        "H_out": c_H_ext,
        "Na_in": c_Na_int,
        "Na_out": c_Na_ext,
    }
    df = pd.DataFrame.from_dict(data_dict, orient="index")
    df.to_csv(filename)


def get_antiporter_df(filename="antiporter_data.csv"):
    df = pd.read_csv(filename, index_col=0)
    return df


def generate_op_fluxes_in_terms_of_tfluxes(G2, G2_df):
    # generate the state probability expressions
    G2_prob_strs = calculations.calc_state_probs(G2, key='name', output_strings=True)
    p1, p2, p3, p4, p5, p6 = G2_prob_strs
    rate_names, sub_map = get_G2_substitution_dict()

    # unfortunately we have to create these here as well
    # so we can use them to read from sub_map
    k12, k21, k23, k32, k34, k43, k45, k54, k56, k65, k61, k16, k14, k41 = sympy.symbols(
        "k12 k21 k23 k32 k34 k43 k45 k54 k56 k65 k61 k16 k14, k41")

    # J_H = j12 - j21
    # J_D = j61 - j16
    J_H_trans = (p1.subs(sub_map) * sub_map[k12] - p2.subs(sub_map) * sub_map[k21]).simplify()
    J_Na_trans = (p6.subs(sub_map) * sub_map[k61] - p1.subs(sub_map) * sub_map[k16]).simplify()

    # convert to lambda functions
    J_H_trans_lambda, J_Na_trans_lambda = expressions.construct_lambda_funcs([J_H_trans, J_Na_trans], rate_names)

    k_leak_vals = [1e-2, 1, 1e2]
    J_H_tflux_arr = np.zeros((len(k_leak_vals),), dtype=np.float64)
    J_Na_tflux_arr = np.zeros((len(k_leak_vals),), dtype=np.float64)
    for i, k_leak in enumerate(k_leak_vals):
        J_H_val = J_H_trans_lambda(
                A_on=G2_df.loc["H_on"].values,
                A_off=G2_df.loc["H_off"].values,
                B_on=G2_df.loc["Na_on"].values,
                B_off=G2_df.loc["Na_off"].values,
                k_AA_H=G2_df.loc["k_AA_H"].values,
                k_AA_Na=G2_df.loc["k_AA_Na"].values,
                k_leak=k_leak,
                A_in=G2_df.loc["H_in"].values,
                A_out=G2_df.loc["H_out"].values,
                B_in=G2_df.loc["Na_in"].values,
                B_out=G2_df.loc["Na_out"].values,
            )
        J_Na_val = J_Na_trans_lambda(
                A_on=G2_df.loc["H_on"].values,
                A_off=G2_df.loc["H_off"].values,
                B_on=G2_df.loc["Na_on"].values,
                B_off=G2_df.loc["Na_off"].values,
                k_AA_H=G2_df.loc["k_AA_H"].values,
                k_AA_Na=G2_df.loc["k_AA_Na"].values,
                k_leak=k_leak,
                A_in=G2_df.loc["H_in"].values,
                A_out=G2_df.loc["H_out"].values,
                B_in=G2_df.loc["Na_in"].values,
                B_out=G2_df.loc["Na_out"].values,
            )
        J_H_tflux_arr[i] = J_H_val[0]
        J_Na_tflux_arr[i] = J_Na_val[0]

    data_arr = np.column_stack((k_leak_vals, J_H_tflux_arr, J_Na_tflux_arr))
    df = pd.DataFrame(data_arr, columns=["k_leak", "J_12", "J_61"])
    return df

if __name__ == "__main__":

    kvals1 = get_G1_rate_matrix()
    kvals2 = get_G2_rate_matrix()
    G1 = generate_input_diagram(kvals1)
    G2 = generate_input_diagram(kvals2)

    cwd = os.getcwd()
    pos = get_node_positions()

    pars = diagrams.generate_partial_diagrams(G2)
    dirpars = diagrams.generate_directional_diagrams(G2)

    G2_cycles = graph_utils.find_all_unique_cycles(G2)
    # change the order of the cycle for the cycles figure (B - C - A)
    G2_cycles_reordered = [G2_cycles[2], G2_cycles[1], G2_cycles[0]]

    G1_net_cycle_flux_func = get_G1_net_cycle_flux_lambda_function(G1)
    G2_net_cycle_flux_funcs = get_G2_net_cycle_flux_lambda_function(G2, G2_cycles=G2_cycles_reordered)

    # create the antiporter dataframe
    create_antiporter_df()
    df = get_antiporter_df()

    leak_arr = np.logspace(-2, 3, 101)

    G1_flux = G1_net_cycle_flux_func(
        A_on=df.loc["H_on"].values,
        A_off=df.loc["H_off"].values,
        B_on=df.loc["Na_on"].values,
        B_off=df.loc["Na_off"].values,
        k_AA_H=df.loc["k_AA_H"].values,
        k_AA_Na=df.loc["k_AA_Na"].values,
        A_in=df.loc["H_in"].values,
        A_out=df.loc["H_out"].values,
        B_in=df.loc["Na_in"].values,
        B_out=df.loc["Na_out"].values,
    )
    G1_flux_arr = np.ones((leak_arr.size,)) * G1_flux

    G2_flux_arr = np.zeros((leak_arr.size, len(G2_net_cycle_flux_funcs)))
    for i, k_leak in enumerate(leak_arr):
        for j, func in enumerate(G2_net_cycle_flux_funcs):
            J_x_G2 = func(
                A_on=df.loc["H_on"].values,
                A_off=df.loc["H_off"].values,
                B_on=df.loc["Na_on"].values,
                B_off=df.loc["Na_off"].values,
                k_AA_H=df.loc["k_AA_H"].values,
                k_AA_Na=df.loc["k_AA_Na"].values,
                k_leak=k_leak,
                A_in=df.loc["H_in"].values,
                A_out=df.loc["H_out"].values,
                B_in=df.loc["Na_in"].values,
                B_out=df.loc["Na_out"].values,
            )
            G2_flux_arr[i, j] = J_x_G2[0]

    # cycle order is B, C, A
    G2_cycle_B_idx = 0
    G2_cycle_prod_idx = 1
    G2_cycle_A_idx = 2
    # assign arrays more useful variable names
    G2_productive_cycle_flux_arr = G2_flux_arr[:, G2_cycle_prod_idx]
    G2_cycle_A_flux_arr = G2_flux_arr[:, G2_cycle_A_idx]
    G2_cycle_B_flux_arr = G2_flux_arr[:, G2_cycle_B_idx]
    # calculate the operational fluxes for cycles A and B for G2
    A_flux_arr = G2_productive_cycle_flux_arr + G2_cycle_A_flux_arr
    B_flux_arr = G2_productive_cycle_flux_arr + G2_cycle_B_flux_arr

    # output plot data to .csv
    plot_data_arr = np.column_stack((leak_arr, G1_flux_arr, A_flux_arr, B_flux_arr))
    plot_df = pd.DataFrame(plot_data_arr, columns=["k_leak", "J (G1)", "J_H+", "J_Na+"])
    plot_df.to_csv("./6wl_analysis_data.csv", index=False)

    # generate the kinetic diagrams for both models
    for diagram, name in zip([G1, G2], ["G", "G_leak"]):
        fname = f"{name}_kinetic_diagram"
        _fig = plotting.draw_diagrams(diagram, pos=pos, font_size=12, curved_arrows=True)
        for file_ext in (".png", ".pdf"):
            _fig.savefig(os.path.join(cwd, fname + file_ext), dpi=300)

    # generate the stacked operational flux figure
    stacked_fig = plot_stacked_flux_fig(leak_arr, G1_flux_arr, G2_productive_cycle_flux_arr, A_flux_arr, B_flux_arr)
    for file_ext in (".png", ".pdf"):
        stacked_fig.savefig(os.path.join(cwd, "6wl_analysis" + file_ext), dpi=300)

    # generate the transition flux data
    tflux_df = generate_op_fluxes_in_terms_of_tfluxes(G2, G2_df=df)
    tflux_df.to_csv("./tflux_data.csv", index=False)

    # other plots that are currently not used
    # plotting.draw_diagrams(pars, pos=pos, panel=True, panel_scale=1.75, rows=1, font_size=12, path=cwd, label="G2_partial_diagrams")
    # plotting.draw_diagrams(dirpars, pos=pos, panel=True, panel_scale=1.75, rows=6, font_size=12, cbt=True, path=cwd, label="G2_directional_diagrams")
    # plotting.draw_cycles(G2, G2_cycles_reordered, pos=pos, panel=True, panel_scale=1.75, font_size=12, path=cwd, label="G2_cycles", curved_arrows=True, cbt=False)
    # cycle_fig = draw_cycles_custom(G2, G2_cycles_reordered, pos=pos, panel_scale=1.75, font_size=12, curved_arrows=True)
