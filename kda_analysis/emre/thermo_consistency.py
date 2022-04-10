"""
Module for verifying the thermodynamic consistency of the sets of
rates used in the analysis of EmrE.
"""

import os
import sys
import numpy as np
import pandas as pd
import networkx as nx

import fig_utils
import operational_flux

from kda import calculations, diagrams, graph_utils


class HiddenPrints:
    """
    Class used to suppress unwanted print statements.
    """

    def __enter__(self):
        self._original_stdout = sys.stdout
        sys.stdout = open(os.devnull, "w")

    def __exit__(self, exc_type, exc_val, exc_tb):
        sys.stdout.close()
        sys.stdout = self._original_stdout


def check_rate_products(K, unique_cycles):
    for cycle in unique_cycles:
        # concatenate a copy of the cycle to the original cycle
        # so pairs of indices loop through the cycle
        extended_cycle = 2 * cycle

        rate_idx = np.zeros((len(cycle), 2), dtype=int)
        for i in range(len(cycle)):
            rate_idx[i] = extended_cycle[i : i + 2]

        forward_rate_product = 1
        reverse_rate_product = 1
        for idx in rate_idx:
            (i, j) = idx
            # the "reverse" rate here is just the opposite direction as the
            # "forward" rate, but it doesn't matter which direction is marked
            # "forward" since we are just comparing the rate products
            forward_rate_product *= K[i, j]
            reverse_rate_product *= K[j, i]

        if not np.isclose(
            forward_rate_product, reverse_rate_product, atol=1e-11, rtol=1e-14
        ):
            raise ValueError(
                f"Forward/reverse rate products do not agree for cycle: {cycle} \n"
                f"Forward rate product: {forward_rate_product} \n"
                f"Reverse rate product: {reverse_rate_product} \n"
            )
    print("--> Passed rate product check")


def check_net_cycle_fluxes(G, dirpar_edges, unique_cycles):
    """
    Verifies that for every cycle in the generated diagram, the calculated net
    cycle flux is very close to zero. Since MultiBind returns a set of rates
    that are at equilibrium, there should be a zero net cycle flux for every
    cycle.
    """
    with HiddenPrints():
        # in order to speed things up we will manually calculate the net
        # cycle fluxes instead of using the built-in KDA function
        # calculations.calc_net_cycle_flux() (so we don't have to calculate
        # sigma every time)
        sigma = calculations.calc_sigma(G, dirpar_edges, key="val")

        for cycle in unique_cycles:
            # normally one has to manually determine the order of the nodes that
            # is considered the "positive" flux direction for a given physiological
            # function to make sure the net cycle flux calculation yields the
            # correct sign, but since we aren't using the net cycle fluxes
            # we can just assign an arbitrary direction and not worry about signs
            order = cycle[:2]
            flux_diags = diagrams.generate_flux_diagrams(G, cycle)
            pi_diff = calculations.calc_pi_difference(
                G,
                cycle,
                order,
                key="val",
            )
            sigma_K = calculations.calc_sigma_K(G, cycle, flux_diags, key="val")
            net_cycle_flux = pi_diff * sigma_K / sigma

            if np.isclose(net_cycle_flux, 0, atol=1e-11, rtol=1e-14):
                continue
            else:
                raise ValueError(
                    f"Calculated a non-zero net cycle flux for cycle: {cycle} \n"
                    f"Calculated net cycle flux: {net_cycle_flux}"
                )
    print("--> Passed net cycle flux check")


def check_thermo_forces(G, cos_dict, name=None):
    for cycle, (order, sign) in cos_dict.items():
        thermo_force = calculations.calc_thermo_force(
            G, cycle, order, key="val", output_strings=False
        )
        if not np.isclose(thermo_force, 0, atol=1e-11, rtol=1e-14):
            raise ValueError(
                f"Calculated thermodynamic force is not zero for cycle {cycle} \n"
                f"Calculated thermodynamic force: {thermo_force}"
            )
    print(f"--> Passed thermodynamic driving force check for {name}")


def check_transition_fluxes(G, prob_arr):
    """
    To the same end as `check_net_cycle_fluxes()`, this verifies the net
    transition fluxes for each pair of states is very close to zero.
    This is done by finding all the unique connections in the diagram G,
    calculating the transition flux in each direction (J_ij and J_ji), then
    checking that they are equal. If they are all equal, then
    J_ij - J_ji = 0 for all i, j, which means we have a set of equilibrium rates.
    """
    unique_edges = diagrams._find_unique_edges(G)
    for edge in unique_edges:
        # assign the tuple values to i and j for readability
        i = edge[0]
        j = edge[1]
        # get the values for the edges from G
        kij = G[i][j][0]["val"]
        kji = G[j][i][0]["val"]
        # get the state probabilities for nodes i and j
        pi = prob_arr[i]
        pj = prob_arr[j]
        # calculate the transition fluxes in both directions
        J_ij = pi * kij
        J_ji = pj * kji
        if not np.isclose(J_ij, J_ji, atol=1e-11, rtol=1e-14):
            raise ValueError(
                f"Forward/reverse transition fluxes do not agree for edge: {edge} \n"
                f"J_({i}, {j}): {J_ij} \n"
                f"J_({j}, {i}): {J_ji} \n"
            )

    print("--> Passed transition fluxes check")


def check_thermo_consistency(G, K, H_cos_dict, D_cos_dict):
    """
    Want to perform a series of checks:
    1.) Check the forward/reverse rate products for each cycle
    2.) Check the net cycle fluxes (with equal concentration ratios) are
        zero for each cycle
    3.) Check that the thermodynamic driving force
        is equal to zero for each cycle
    4.) if we have the state probabilities, we can also check the transition
        fluxes for each edge pair, but probably not necessary
    """
    unique_cycles = graph_utils.find_all_unique_cycles(G)
    dirpar_edges = diagrams.generate_directional_partial_diagrams(G, return_edges=True)
    state_probs = calculations.calc_state_probs_from_diags(
        G=G, dirpar_edges=dirpar_edges, key="val"
    )

    check_rate_products(K=K, unique_cycles=unique_cycles)
    check_net_cycle_fluxes(G=G, dirpar_edges=dirpar_edges, unique_cycles=unique_cycles)
    check_thermo_forces(G=G, cos_dict=H_cos_dict, name="H+")
    check_thermo_forces(G=G, cos_dict=D_cos_dict, name="Substrate")
    check_transition_fluxes(G=G, prob_arr=state_probs)


def run_thermo_consistency_checks(df, data_path=None):
    H_edges = [(0, 6, 0), (1, 7, 0)]
    D_edges = [(3, 5, 0), (1, 7, 0)]

    n_datasets = df.shape[0]
    is_valid_dataset = np.ones(n_datasets, dtype=bool)

    for i, row in enumerate(df.iterrows()):
        fig_key = row[1]["figure"]
        k_data = row[1].values

        K = fig_utils.get_K(
            k_data[0],
            k_data[1],
            k_data[2],
            k_data[3],
            k_data[4],
            k_data[5],
            k_data[6],
            k_data[7],
            k_data[8],
            k_data[9],
            k_data[10],
            k_data[11],
            k_data[12],
            k_data[13],
            k_data[14],
            k_data[15],
            k_data[16],
            k_data[17],
            k_data[18],
            k_data[19],
            k_data[20],
            k_data[21],
            k_data[22],
            k_data[23],
        )

        G = nx.MultiDiGraph()
        graph_utils.generate_edges(G, K)

        H_cycles_valid, H_cycles_invalid = operational_flux.find_relevant_cycles(
            G, H_edges
        )
        D_cycles_valid, D_cycles_invalid = operational_flux.find_relevant_cycles(
            G, D_edges
        )
        H_cos_dict, D_cos_dict = operational_flux.get_cos_dicts(
            H_cycles_valid, D_cycles_valid
        )

        print(f"--> Running thermo consistency check {i+1}/{n_datasets}")
        sys.stdout.flush()
        try:
            check_thermo_consistency(
                G=G, K=K, H_cos_dict=H_cos_dict, D_cos_dict=D_cos_dict
            )
        except ValueError as err:
            print(err)
            # if it does not pass, mark invalid
            is_valid_dataset[i] = False
            continue

    # add column for thermodynamic consistency
    n_valid_datasets = np.count_nonzero(is_valid_dataset)
    print(f"\n--> Valid datasets: {n_valid_datasets}/{n_datasets}")
    if not data_path is None:
        df["thermodynamically consistent"] = is_valid_dataset
        df.to_csv(data_path, index=False)


if __name__ == "__main__":
    data_path = "data/figure_data.csv"
    df = pd.read_csv(data_path)
    run_thermo_consistency_checks(df=df, data_path=data_path)
