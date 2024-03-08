import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import tflux_utils


def calc_proton_operational_flux(prob_arr, k_AA_anti, k_AA_sym, option):
    """
    Proton Flux Equivalencies:
    J_H  =  J_{7,1} + J_{8,2}  =  J_{1,3} + J_{2,4}  =  J_{6,8} + J_{5,7}

    Net Transition Flux Definition:
    J_{i,j} = (k_{i,j} * p_i) - (k_{j,i} * p_j)

    Rate definitions:
        k31: H_on * H_out
        k13: H_off
        k57: H_on * H_in
        k75: H_off
        k42: H_on * H_out
        k24: H_off
        k68: H_on * H_in
        k86: H_off
        k17: k_AA_anti
        k71: k_AA_anti
        k28: k_AA_sym
        k82: k_AA_sym

    """
    H_on = 1e10
    H_off = 1e3
    H_in = 10 ** (-6.5)
    H_out = 10 ** (-7.5)

    # calculate the transition fluxes
    if option == "A":
        # J_71 = (k71 * p_7) - (k17 * p_1)
        J_1 =  (k_AA_anti * prob_arr[6]) - (k_AA_anti * prob_arr[0])
        # J_82 = (k82 * p_8) - (k28 * p_2)
        J_2 =  (k_AA_sym * prob_arr[7]) - (k_AA_sym * prob_arr[1])
    elif option == "B":
        # J_13 = (k13 * p_1) - (k31 * p_3)
        J_1 = (H_off * prob_arr[0]) - (H_on * H_out * prob_arr[2])
        # J_24 = (k24 * p_2) - (k42 * p_4)
        J_2 = (H_off * prob_arr[1]) - (H_on * H_out * prob_arr[3])
    elif option == "C":
        # J_68 = (k68 * p_6) - (k86 * p_8)
        J_1 =  (H_on * H_in * prob_arr[5]) - (H_off * prob_arr[7])
        # J_57 = (k57 * p_5) - (k75 * p_7)
        J_2 =  (H_on * H_in * prob_arr[4]) - (H_off * prob_arr[6])
    else:
        raise ValueError("Invalid option, please choose A, B, or C.")


    # calculate the operational fluxes
    J_H_arr = J_1 + J_2

    return J_H_arr


def calc_drug_operational_flux(prob_arr, k_AA_anti, k_AA_sym, option):
    """
    Drug Flux Equivalencies:
    J_D  =  J_{6,4} + J_{8,2}  =  J_{7,8} + J_{5,6}  =  J_{2,1} + J_{4,3}

    Net Transition Flux Definition:
    J_{i,j} = (k_{i,j} * p_i) - (k_{j,i} * p_j)

    Rate Definitions:
        k34: D_on * c_D
        k43: D_off
        k56: D_on * c_D
        k65: D_off
        k12: D_on * c_D
        k21: D_off
        k78: D_on * c_D
        k87: D_off
        k46: k_AA_anti
        k64: k_AA_anti
        k28: k_AA_sym
        k82: k_AA_sym
    """
    D_on = 1e7
    D_off = 1
    c_D = 25e-9

    # calculate the transition fluxes
    if option == "A":
        # J_64 = (k64 * p_6) - (k46 * p_4)
        J_1 = (k_AA_anti * prob_arr[5]) - (k_AA_anti * prob_arr[3])
        # J_82 = (k82 * p_8) - (k28 * p_2)
        J_2 = (k_AA_sym * prob_arr[7]) - (k_AA_sym * prob_arr[1])
    elif option == "B":
        # J_78 = (k78 * p_7) - (k87 * p_8)
        J_1 = (D_on * c_D * prob_arr[6]) - (D_off * prob_arr[7])
        # J_56 = (k56 * p_5) - (k65 * p_6)
        J_2 = (D_on * c_D * prob_arr[4]) - (D_off * prob_arr[5])
    elif option == "C":
        # J_21 = (k21 * p_2) - (k12 * p_1)
        J_1 = (D_off * prob_arr[1]) - (D_on * c_D * prob_arr[0])
        # J_43 = (k43 * p_4) - (k34 * p_3)
        J_2 = (D_off * prob_arr[3]) - (D_on * c_D * prob_arr[2])
    else:
        raise ValueError("Invalid option, please choose A, B, or C.")

    # calculate the operational fluxes
    J_D_arr = J_1 + J_2

    return J_D_arr


def plot_probabilities(RAA_arr, prob_arr, fig_key):
    fig = plt.figure(figsize=(3.25, 4.25), tight_layout=True)
    ax = fig.add_subplot(211)

    for i in range(prob_arr.shape[0]):
        _label = r"$p_{%d}$" % (i + 1)
        ax.semilogx(
            RAA_arr,
            prob_arr[i],
            ls="-",
            label=_label,
        )

    ax.axvline(
        x=1,
        ymin=0,
        ymax=1,
        ls="--",
        color="black",
    )
    ax.legend(loc="upper left", bbox_to_anchor=(1, 1))

    ax.set_title(r"EmrE State Probabilities vs. $R_${AA}")
    ax.set_xlabel(r"$R_{AA}$")
    ax.set_ylabel("Probability")


def plot_operational_fluxes(RAA_arr, k_AA_anti, k_AA_sym, prob_arr, fig_key, option):
    """
    Want to plot the operational fluxes for H & D for EmrE
    using the transition fluxes.

        J_H  =  J_{7,1} + J_{8,2}  =  J_{1,3} + J_{2,4}  =  J_{6,8} + J_{5,7} 
        J_D  =  J_{6,4} + J_{8,2}  =  J_{7,8} + J_{5,6}  =  J_{2,1} + J_{4,3}

    """

    # calculate the operational fluxes
    J_H_arr = calc_proton_operational_flux(prob_arr, k_AA_anti=k_AA_anti, k_AA_sym=k_AA_sym, option=option)
    J_D_arr = calc_drug_operational_flux(prob_arr, k_AA_anti=k_AA_anti, k_AA_sym=k_AA_sym, option=option)

    ######################################################
    ### Operational Flux Figure ##########################
    ######################################################

    fig = plt.figure(figsize=(6, 8), tight_layout=True)
    ax = fig.add_subplot(211)
    ax_stoich = fig.add_subplot(212)

    for _ax in (ax, ax_stoich):
        _ax.grid(True)

    H_colour = "#41b6c4"
    D_colour = "#a1dab4"

    ax.semilogx(
        RAA_arr,
        J_H_arr,
        ls="-",
        color=H_colour,
        label=r"$J_{\mathrm{H}^{+}}$",
    )

    ax.semilogx(
        RAA_arr,
        J_D_arr,
        ls="-",
        color=D_colour,
        label=r"$J_{\mathrm{D}}$",
    )

    ax.axvline(
        x=1,
        ymin=0,
        ymax=1,
        ls="--",
        color="black",
    )
    ax.legend(loc="upper left", bbox_to_anchor=(1, 1))

    ######################################################
    ### Stoichiometry Figure #############################
    ######################################################

    ax_stoich.semilogx(
        RAA_arr,
        J_D_arr/J_H_arr,
        ls="-",
        color="grey",
        label=r"$J_{\mathrm{H}^{+}}$",
    )

    ax_stoich.axvline(
        x=1,
        ymin=0,
        ymax=1,
        ls="--",
        color="black",
    )

    ax.set_ylabel("Flux (" + r"$\mathrm{s}^{-1}$" + ")")
    ax_stoich.set_xlabel(r"$R_{AA}$")


def plot_transition_flux(R, edges, fluxes):

    fig = plt.figure(figsize=(4, 3), tight_layout=True)
    ax = fig.add_subplot(111)

    color_list = [
        "#e6194B", "#3cb44b", "#ffe119", "#4363d8", 
        "#f58231", "#911eb4", "#42d4f4", "#f032e6", 
        "#bfef45", "#fabed4", "#469990", "#9A6324", 
        "#800000", "#808000", "#000075", "#a9a9a9", 
        "#000000",
    ]
    color_idx = 0

    def get_linestyle(idx):
        if idx % 2 == 0:
            return "-"
        else:
            return "--"

    for edge, flux in zip(edges, fluxes):
        label = r"$J_{%d, %d}$" % (edge[0]+1,edge[1]+1)
        ax.semilogx(
            R, 
            flux, 
            color=color_list[color_idx], 
            # ls=get_linestyle(color_idx),
            lw=1,
            label=label,
        )
        color_idx += 1
    ax.set_title("EmrE Net Transition Fluxes")
    ax.set_ylabel(r"Flux (s$^{-1}$)")
    ax.set_xlabel(r"$R_{\mathrm{AA}}$")

    ax.axvline(x=1e0, ymin=0, ymax=1, ls="--", color="black")
    plt.legend(bbox_to_anchor=(1, 1), loc="upper left", ncol=1)
    ax.grid(True)

    plt.close()
    return fig


def main():
   # collect the rate matrix and generate a MultiDiGraph for EmrE
    K, G = tflux_utils.get_K_and_G()
    expressions = tflux_utils.generate_probability_expressions(G)

    for fig_key in ("7A", "7B"):

        sympy_funcs = tflux_utils.generate_sympy_funcs(expressions, fig_key=fig_key)
        pi_lambda = tflux_utils.generate_lambda_funcs(sympy_funcs, fig_key=fig_key)

        # plug in values and create plot arrays

        if fig_key == "7A":
            # don't round values since some of the array elements end up
            # being the same value, and generate 101 elements so we hit the
            # integers 1, 10, and 100 precisely
            n_datasets = 201
            # favor antiport initially
            k_AA_anti = np.logspace(-3, 5, n_datasets)[::-1]
            k_AA_sym = np.logspace(-3, 5, n_datasets)

            RAA_arr = k_AA_anti/k_AA_sym

            # create an array to store all state probability values
            # for `n_datasets`
            prob_arr = np.zeros((G.number_of_nodes(), n_datasets), dtype=np.float64)
            for i, p_i in enumerate(pi_lambda):
                for j, (k_anti, k_sym) in enumerate(zip(k_AA_anti, k_AA_sym)):
                    prob_arr[i, j] = p_i(k_AA_anti=k_anti, k_AA_sym=k_sym)

            # create raw data for operational flux (in terms of net transition fluxes)
            # table in the KDA paper. We have to use option "A" here since we are using
            # J_H = J_{7,1} + J_{8,2} and J_D = J_{6,4} + J_{8,2} in the paper
            J_H_arr = calc_proton_operational_flux(prob_arr, k_AA_anti=k_AA_anti, k_AA_sym=k_AA_sym, option="A")
            J_D_arr = calc_drug_operational_flux(prob_arr, k_AA_anti=k_AA_anti, k_AA_sym=k_AA_sym, option="A")
            df = pd.DataFrame(
                np.column_stack((RAA_arr, J_H_arr, J_D_arr)),
                columns=["R_AA", "J_H = J_{7,1} + J_{8,2}", "J_D = J_{6,4} + J_{8,2}"],
                )

            # code to generate the net transition fluxes
            # as a function of RAA. Commented out since it was
            # never used in the paper.
            # edges = []
            # for i, j, k in G.edges:
            #     if (j, i) not in edges:
            #         edges.append((i, j))

            # tflux_arr = np.zeros((len(edges), n_datasets), dtype=np.float64)
            # for e, edge in enumerate(edges):
            #     i, j = edge
            #     for v, (k_anti, k_sym) in enumerate(zip(k_AA_anti, k_AA_sym)):
            #         K_tmp = tflux_utils.get_K_vals(k_AA_anti=k_anti, k_AA_sym=k_sym)
            #         # J_{i,j} = (k_{i,j} * p_i) - (k_{j,i} * p_j)
            #         j_ij = K_tmp[i, j] * pi_lambda[i](k_anti, k_sym)
            #         j_ji = K_tmp[j, i] * pi_lambda[j](k_anti, k_sym)
            #         tflux_arr[e, v] = j_ij - j_ji

            # fig = plot_transition_flux(R=RAA_arr, edges=edges, fluxes=tflux_arr)
            # for ext in ("png", "pdf"):
            #     fig.savefig(f"transition_fluxes_{fig_key}.{ext}", dpi=300)

        elif fig_key == "7B":
            n_datasets = 201

            # favor antiport initially by using large EHD values first
            k_EH_H_arr = np.logspace(-2, 8, n_datasets)
            k_EHD_H_arr = np.logspace(-2, 8, n_datasets)[::-1]
            k_ED_D_arr = np.logspace(-4, 6, n_datasets)
            k_EHD_D_arr = np.logspace(-4, 6, n_datasets)[::-1]
            # only calculate for a single case of k_AA for simplicity
            # k_AA = np.array([1, 10, 100, 1000])
            k_AA = 100

            R_off_arr = k_EH_H_arr / k_EHD_H_arr

            # create an array to store all state probability values
            # for `n_datasets`
            prob_arr = np.zeros((G.number_of_nodes(), n_datasets), dtype=np.float64)
            for i, p_i in enumerate(pi_lambda):
                for j, (k_EH_H, k_EHD_H, k_ED_D, k_EHD_D) in enumerate(zip(k_EH_H_arr, k_EHD_H_arr, k_ED_D_arr, k_EHD_D_arr)):
                    prob_arr[i, j] = p_i(k_EH_H=k_EH_H, k_EHD_H=k_EHD_H, k_ED_D=k_ED_D, k_EHD_D=k_EHD_D, k_AA=k_AA)

            # J_71 = (k71 * p_7) - (k17 * p_1)
            J_71 =  (k_AA * prob_arr[6]) - (k_AA * prob_arr[0])
            # J_82 = (k82 * p_8) - (k28 * p_2)
            J_82 =  (k_AA * prob_arr[7]) - (k_AA * prob_arr[1])
            # J_64 = (k64 * p_6) - (k46 * p_4)
            J_64 = (k_AA * prob_arr[5]) - (k_AA * prob_arr[3])

            # create raw data for operational flux (in terms of net transition fluxes)
            # table in the KDA paper. We have to use option "A" here since we are using
            # J_H = J_{7,1} + J_{8,2} and J_D = J_{6,4} + J_{8,2} in the paper
            J_H_arr = J_71 + J_82
            J_D_arr = J_64 + J_82

            df = pd.DataFrame(
                np.column_stack((R_off_arr, J_H_arr, J_D_arr)),
                columns=["R_off", "J_H = J_{7,1} + J_{8,2}", "J_D = J_{6,4} + J_{8,2}"]
            )

        # output dataframe as a .csv
        df.to_csv(f"./opflux_from_tflux_{fig_key}.csv", index=False)

if __name__ == "__main__":
    main()
