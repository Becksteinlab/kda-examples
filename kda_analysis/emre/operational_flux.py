import os
import sys

# have to raise the recursion limit or SymPy raises a RecursionError when
# attempting to parse the EmrE 8-state model expression for Sigma
sys.setrecursionlimit(5000)
import time

import dill
import numpy as np
from sympy import simplify, lambdify
from sympy.parsing.sympy_parser import parse_expr
from kda import calculations, diagrams, graph_utils


def find_relevant_cycles(G, edges):
    """
    Finds the cycles in the input diagram G that contain only 1 edge in the
    input list of edges. This is used to determine which cycles contribute to
    the operational flux of protons/drug molecules in the 8-state model of EmrE.

    Parameters
    ----------
    G : NetworkX MultiDiGraph Object
        Diagram for the 8-state model of EmrE.
    edges : list of edge 3-tuples
        List of edge tuples, where each edge is a contributing process to the
        flux of a given species (proton/drug) For example, in the EmrE
        8-state model, the 2 edges that correspond to the transport of drug
        molecules are `[(1, 7, 0), (3, 5, 0)]`. Note: these have to be 3-tuples
        or the comparison operations will not work.

    Returns
    -------
    (valid_cycles, invalid_cycles) : tuple of nested lists
        A tuple of the lists of cycles that contribute to the operational
        flux (valid) as well as the non-contributing cycles (invalid). Each
        cycle in each cycle list is a list containing the integers for the
        nodes in that cycle.
    """
    # get all cycles in the input diagram
    cycle_list = graph_utils.find_all_unique_cycles(G)

    # create two lists, one for operational flux contributors (valid)
    # cycles and another for non-contributors (invalid)
    valid_cycles = []
    invalid_cycles = []

    for cycle in cycle_list:
        # get the unidirectional/unique cycle edges
        unidirectional_edges = diagrams._construct_cycle_edges(cycle)
        # use the inidirectional edges to obtain the bidirectional edge list
        all_cycle_edges = diagrams._append_reverse_edges(unidirectional_edges)

        # count how many input edges are in the given cycle
        count = 0
        for edge in edges:
            if edge in all_cycle_edges:
                count += 1
        # if there is exactly 1 edge in the cycle, this cycle
        # is a contributor to the operational flux
        if count == 1:
            if not cycle in valid_cycles:
                # only add to valid cycles if it is not already in valid cycles
                valid_cycles.append(cycle)
        elif count == 0:
            # if there are zero edges in the bidirectional edge list
            # this is not a contributing cycle
            if not cycle in invalid_cycles:
                invalid_cycles.append(cycle)
        elif count == 2:
            # if there are 2 edges in the bidirectional edge list, this can mean
            # many things. For the EmrE 8-state model, all of these cases happen
            # to be in a configuration where a given species
            # (proton/drug molecule) is being moved in and out in the same
            # cycle, thus not contributing to the operational flux. For more
            # complex models this may not be the case.
            if not cycle in invalid_cycles:
                invalid_cycles.append(cycle)
        elif count > 2:
            # For this particular case (EmrE 8-state model) these cases
            # do not show up, so if they are observed it is due to an error
            raise NotImplementedError(
                f"Counted {count} edges in cycle {cycle} from \n"
                f"input cycle edges: {edges}"
            )

    # check that the cycles are conserved
    total_cycles = len(cycle_list)
    counted_cycles = len(valid_cycles) + len(invalid_cycles)
    if counted_cycles != total_cycles:
        raise ValueError(
            f"Failed to sort all cycles. \n"
            f"{counted_cycles} accounted for out of {total_cycles}"
        )

    return valid_cycles, invalid_cycles


def get_cos_dicts(H_cycles_valid, D_cycles_valid):
    """
    Generates the Proton and Drug Cycle/Order/Sign dictionaries for the
    8-state EmrE model.

    Parameters
    ----------
    H_cycles_valid : nested list
        List of all cycles (in the 8-state model) that have a net transport
        of 1 proton per cycle completion. The cycles are generated using
        `find_relevant_cycles()`.

    D_cycles_valid : nested list
        List of all cycles (in the 8-state model) that have a net transport
        of 1 drug molecule per cycle completion. The cycles are generated using
        `find_relevant_cycles()`.

    Returns
    -------
    (H_cos_dict, D_cos_dict) : tuple of dictionaries
        Tuple of Cycle/Order/Sign dictionaries.

    Examples
    --------
    The fully assembled Cycle/Order/Sign dictionary for the operational flux
    of protons in the EmrE 8-state model:

    cos_dict = {
        (0, 6, 7, 5, 4, 2, 3, 1): ([6, 0], '+'),
        (0, 6, 7, 5, 4, 2): ([6, 0], '+'),
        (0, 6, 7, 5, 3, 2): ([6, 0], '+'),
        (0, 6, 7, 5, 3, 1): ([6, 0], '+'),
        (0, 6, 4, 5, 3, 2): ([6, 0], '+'),
        (0, 6, 4, 5, 3, 1): ([6, 0], '+'),
        (0, 6, 4, 2, 3, 1): ([6, 0], '+'),
        (0, 6, 4, 2): ([6, 0], '+'),
        (0, 2, 4, 6, 7, 1): ([0, 2], '+'),
        (0, 2, 4, 5, 7, 1): ([0, 2], '+'),
        (0, 2, 3, 5, 7, 1): ([0, 2], '+'),
        (0, 2, 3, 5, 4, 6, 7, 1): ([0, 2], '+'),
        (1, 7, 6, 4, 5, 3): ([1, 3], '+'),
        (1, 7, 6, 4, 2, 3): ([1, 3], '+'),
        (1, 7, 5, 4, 2, 3): ([1, 3], '+'),
        (1, 7, 5, 3): ([1, 3], '+'),
    }
    """
    # create a list of minus/plus sign strings
    minus = ["-"]
    plus = ["+"]
    # all of the proton cycles are positive contributors
    H_signs = 16 * plus
    # the first half of the drug cycles are negative
    # contributors, the rest are positive contributors
    D_signs = 8 * minus + 8 * plus
    # each cycle order is a list of nodes (a subset of the cycle itself)
    # that describe which direction is CCW. These were determined via
    # visual inspection for each case.
    H_cycle_order = [
        [6, 0],
        [6, 0],
        [6, 0],
        [6, 0],
        [6, 0],
        [6, 0],
        [6, 0],
        [6, 0],
        [0, 2],
        [0, 2],
        [0, 2],
        [0, 2],
        [1, 3],
        [1, 3],
        [1, 3],
        [1, 3],
    ]
    D_cycle_order = [
        [6, 0],
        [6, 0],
        [6, 0],
        [6, 0],
        [6, 0],
        [6, 0],
        [6, 0],
        [6, 0],
        [2, 4],
        [2, 4],
        [2, 4],
        [2, 4],
        [2, 4],
        [2, 4],
        [2, 4],
        [2, 4],
    ]

    # assemble the Cycle/Order/Sign dictionaries
    H_cos_dict = {}
    for cycle, order, sign in zip(H_cycles_valid, H_cycle_order, H_signs):
        H_cos_dict[tuple(cycle)] = (order, sign)
    D_cos_dict = {}
    for cycle, order, sign in zip(D_cycles_valid, D_cycle_order, D_signs):
        D_cos_dict[tuple(cycle)] = (order, sign)
    return H_cos_dict, D_cos_dict


def get_op_cycle_flux_funcs(
    G,
    cos_dict,
    sub_dict,
    rate_names=None,
    key="name",
):
    """
    Generates the Python lambda function and SymPy expression for the
    operational flux of a species (proton/drug) of the EmrE 8-state model.

    Parameters
    ----------
    G : NetworkX MultiDiGraph Object
        Input diagram.
    cos_dict : dict
        Tuple of Cycle/Order/Sign dictionaries, see `get_cos_dicts()` for
        a full description and example.
    sub_dict : dict
        A dictionary of substitutions to make for SymPy expression manipulation,
        where the keys are the current variable names, the values are the
        new variable names.
    rate_names : list
        List of variable names for the Sympy to Python
        lambda function conversion.
    key : str
        The key for the input diagram used for retrieving the string names of
        the rate matrix variables. Default is "name"

    Returns
    -------
    (lambdify_func, simplified_func) : tuple
        A tuple of the Python lambda function and the SymPy function (simplified
        with the `sub_dict`) for the operational flux, where the lambda
        function was generated from the simplified SymPy expression.
    """
    # get the directional partial diagram (edges)
    dir_par_edges = diagrams.generate_directional_diagrams(G, return_edges=True)
    # use the directional partial diagram info to calculate sigma (normalization)
    sigma_str = calculations.calc_sigma(G, dir_par_edges, key=key, output_strings=True)
    # initialize an empty list for storing net cycle flux numerators
    numerators = []

    # iterate over the input cycles
    for cycle, (order, sign) in cos_dict.items():

        if sign == "+":
            # if the sign is positive, retrieve the Pi difference equation
            # for the cycle
            pi_diff_str = calculations.calc_pi_difference(
                G, cycle, order, key, output_strings=True
            )
        elif sign == "-":
            # if the sign is negative, retrieve the Pi difference equation
            # for the reverse ordered cycle
            rev_order = order[::-1]
            pi_diff_str = calculations.calc_pi_difference(
                G, cycle, rev_order, key, output_strings=True
            )

        # generate the flux diagrams for the cycle
        flux_diags = diagrams.generate_flux_diagrams(G, cycle)
        # use the flux diagrams to calculate `sigma_K`
        sigma_K_str = calculations.calc_sigma_K(
            G, cycle, flux_diags, key, output_strings=True
        )

        if not isinstance(sigma_K_str, str):
            # if 1 is returned for `sigma_k`, our numerator is just
            # the Pi difference
            numerators.append(pi_diff_str)
        else:
            # if `sigma_K` is non-trivial, construct the numerator from
            # both pieces
            numerators.append(str(parse_expr(pi_diff_str) * parse_expr(sigma_K_str)))

    # sum all of the numerators
    num_sum = "+".join(numerators)
    # normalize the numerator with `sigma`
    func = parse_expr(num_sum) / parse_expr(sigma_str)

    # now use the input substitution dictionary to simplify the final expression
    simplified_func = simplify(func.subs(sub_dict))
    if rate_names is None:
        return simplified_func
    else:
        # construct the python/lambda function for calculations
        lambdify_func = lambdify(rate_names, simplified_func, "numpy")
        return lambdify_func, simplified_func


def get_norm_func(G, sub_dict, rate_names, key="name", return_sympy=None):
    """
    Constructs the analytic normalization function (`Sigma`) in SymPy or
    Python lambda function form.

    Parameters
    ----------
    G : NetworkX MultiDiGraph Object
        Input diagram.
    sub_dict : dict
        A dictionary of substitutions to make for SymPy expression manipulation,
        where the keys are the current variable names, the values are the
        new variable names.
    rate_names : list
        List of variable names for the Sympy to Python
        lambda function conversion.
    key : str
        The key for the input diagram used for retrieving the string names of
        the rate matrix variables. Default is "name"
    return_sympy : bool
        Binary used to determine whether to return the Python lambda function
        or the SymPy expression. Default is to return the lambda function.

    Returns
    -------
    norm_lf/norm_sf : Python lambda function/SymPy function
        Python lambda function or SymPy function (simplified with the
        `sub_dict`) for the normalization function, `Sigma`, where the lambda
        function was generated from the simplified SymPy expression.
    """
    start = time.time()
    dir_par_edges = diagrams.generate_directional_partial_diagrams(G, return_edges=True)
    sigma_str = calculations.calc_sigma(G, dir_par_edges, key, output_strings=True)
    norm_f = parse_expr(sigma_str)
    norm_sf = simplify(norm_f.subs(sub_dict))
    if return_sympy:
        end = time.time() - start
        print(f"--> SymPy normalization function generation completed in {end:.2f} s")
        return norm_sf
    else:
        norm_lf = lambdify(rate_names, norm_sf, "numpy")
        end = time.time() - start
        print(f"--> Lambda normalization function generation completed in {end:.2f} s")
        return norm_lf


def get_single_cycle_flux_numerators(
    G,
    cos_dict,
    sub_dict,
    rate_names,
    key="name",
):
    """
    Generates a list of net cycle flux functions for the 8-state model of EmrE.
    These need to be normalized by using the function created in
    `get_norm_func()`, and can be summed to calculate the operational flux
    (when properly normalized).

    Parameters
    ----------
    G : NetworkX MultiDiGraph Object
        Input diagram.
    cos_dict : dict
        Tuple of Cycle/Order/Sign dictionaries, see `get_cos_dicts()` for
        a full description and example.
    sub_dict : dict
        A dictionary of substitutions to make for SymPy expression manipulation,
        where the keys are the current variable names, the values are the
        new variable names.
    rate_names : list
        List of variable names for the Sympy to Python
        lambda function conversion.
    key : str
        The key for the input diagram used for retrieving the string names of
        the rate matrix variables. Default is "name"

    Returns
    -------
    lambda_funcs : list of Python lambda functions
        List of the Python lambda functions for each numerator (net cycle flux)
        that contributes to the operational flux of a species (proton/drug) for
        the EmrE 8-state model.
    """
    lambda_funcs = []
    for cycle, (order, sign) in cos_dict.items():
        cycle = list(cycle)
        if sign == "-":
            # if there is a negative sign, flip the cycle order
            order = order[::-1]
        pi_diff_str = calculations.calc_pi_difference(
            G=G, cycle=cycle, order=order, key=key, output_strings=True
        )
        flux_diags = diagrams.generate_flux_diagrams(G=G, cycle=cycle)
        sigma_K_str = calculations.calc_sigma_K(
            G=G, cycle=cycle, flux_diags=flux_diags, key=key, output_strings=True
        )
        if not isinstance(sigma_K_str, str):
            func = pi_diff_str
        else:
            func = str(parse_expr(pi_diff_str) * parse_expr(sigma_K_str))

        # run the function through the SymPy parser
        func = parse_expr(func)
        # simplify the expression using the substitution dictionary
        simplified_func = simplify(func.subs(sub_dict))
        # convert the SymPy function into a Python lambda function
        lambda_func = lambdify(rate_names, simplified_func, "numpy")
        lambda_funcs.append(lambda_func)

    return lambda_funcs


def get_cycle_order(cycle_idx):
    # manually set the order to CCW
    if cycle_idx < 14:
        order = [6, 0]
    elif (cycle_idx >= 14) and (cycle_idx < 21):
        order = [0, 2]
    elif (cycle_idx >= 21) and (cycle_idx < 25):
        order = [7, 1]
    elif (cycle_idx >= 25) and (cycle_idx < 27):
        order = [2, 4]
    elif cycle_idx == 27:
        order = [4, 6]
    else:
        raise ValueError(f"Too many cycles detected. Expected 28, detected {len(all_cycles)}.")
    return order


def generate_net_cycle_flux_sympy_funcs(G, sub_dict, dir_path):
    all_cycles = graph_utils.find_all_unique_cycles(G)
    dir_edges = diagrams.generate_directional_diagrams(G, return_edges=True)
    sigma_str = calculations.calc_sigma(G, dir_edges, key="name", output_strings=True)

    data_dict = {}
    for i, cycle in enumerate(all_cycles):

        filepath = os.path.join(dir_path, f"fig_7A_cycle_{i+1}_net_cycle_flux.txt")

        if not os.path.isfile(filepath):
            print(
                f"No SymPy function found at location {filepath} \n"
                f"Generating new sympy function..."
            )
            print("=" * 20)
            print(f"Cycle {i+1}: {cycle}")
            order = get_cycle_order(cycle_idx=i)
            print(f"Cycle {i+1} order: {order}")
            sys.stdout.flush()
            flux_diags = diagrams.generate_flux_diagrams(G, cycle=cycle)
            pi_diff_str = calculations.calc_pi_difference(G, cycle, order=order, key="name", output_strings=True)
            sigma_K_str = calculations.calc_sigma_K(G, cycle, flux_diags, key="name", output_strings=True)

            if sigma_K_str == 1:
                numerator = pi_diff_str
            else:
                numerator = str(parse_expr(pi_diff_str) * parse_expr(sigma_K_str))

            sympy_func = parse_expr(numerator) / parse_expr(sigma_str)
            sympy_func = simplify(sympy_func.subs(sub_dict))
            dill.dump(sympy_func, open(filepath, "wb"))


def get_net_cycle_flux_cycles_and_funcs(G, rate_names, dir_path):
    all_cycles = graph_utils.find_all_unique_cycles(G)

    data_dict = {}
    for i, cycle in enumerate(all_cycles):
        filepath = os.path.join(dir_path, f"fig_7A_cycle_{i+1}_net_cycle_flux.txt")
        sympy_func = dill.load(open(filepath, "rb"))
        lambdify_func = lambdify(rate_names, sympy_func, "numpy")

        cycle_key = f"Cycle {i+1}"
        data_dict[cycle_key] = {}
        data_dict[cycle_key]["cycle"] = cycle
        data_dict[cycle_key]["func"] = lambdify_func

    return data_dict

# def get_operational_thermo_force(G, cos_dict, sub_dict):
#     """
#
#     Parameters
#     ----------
#
#
#     Returns
#     -------
#
#     """
#     forces = []
#     for cycle, (order, sign) in cos_dict.items():
#         f = calculations.calc_thermo_force(
#             G, cycle, order, key="name", output_strings=True
#         )
#         f_simp = f.subs(sub_dict)
#         forces.append(f_simp)
#     return forces


# def get_analytic_thermo_force_func(forces, cos_dict):
#     """
#
#     Parameters
#     ----------
#
#
#     Returns
#     -------
#
#     """
#     pos = []
#     neg = []
#     for i, (cycle, (order, sign)) in enumerate(cos_dict.items()):
#         if sign == "+":
#             pos.append(forces[i])
#         elif sign == "-":
#             neg.append(forces[i])
#         else:
#             print("Sign could not be determined.")
#     pos_str = np.sum(pos)
#     neg_str = np.sum(neg)
#     force_str = pos_str - neg_str
#     return force_str


# def get_transition_flux_numerators(
#     G, edge_cycles, cos_dict, sub_dict, rate_names, key="name"
# ):
#     """
#
#     Parameters
#     ----------
#
#
#     Returns
#     -------
#
#     """
#     edge_flux_funcs = []
#     for cycle_list in edge_cycles:
#         cycle_numerators = []
#         for cycle in cycle_list:
#             order, sign = cos_dict[tuple(cycle)]
#             if sign == "+":
#                 pi_diff_str = calculations.calc_pi_difference(
#                     G, cycle, order, key, output_strings=True
#                 )
#             elif sign == "-":
#                 rev_order = order[::-1]
#                 pi_diff_str = calculations.calc_pi_difference(
#                     G, cycle, rev_order, key, output_strings=True
#                 )
#             flux_diags = diagrams.generate_flux_diagrams(G, cycle)
#             sigma_K_str = calculations.calc_sigma_K(
#                 G, cycle, flux_diags, key, output_strings=True
#             )
#             if not isinstance(sigma_K_str, str):
#                 cycle_numerators.append(pi_diff_str)
#             else:
#                 cycle_numerators.append(
#                     str(parse_expr(pi_diff_str) * parse_expr(sigma_K_str))
#                 )
#         num_sum = "+".join(cycle_numerators)
#         func = parse_expr(num_sum)
#         simplified_func = simplify(func.subs(sub_dict))
#         lambdify_func = lambdify(rate_names, simplified_func, "numpy")
#         edge_flux_funcs.append(lambdify_func)
#     return edge_flux_funcs
