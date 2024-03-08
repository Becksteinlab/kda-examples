import os
import sys
# have to raise the recursion limit or SymPy raises a RecursionError when
# attempting to parse the EmrE 8-state model expression for Sigma
sys.setrecursionlimit(5000)

import dill
dill.settings["recursive"] = True

import numpy as np
import networkx as nx
from sympy import symbols, simplify, lambdify
from sympy.parsing.sympy_parser import parse_expr

from kda import graph_utils, calculations


def get_K(
    k31=1, k13=1, k57=1, k75=1, k42=1, k24=1, 
    k68=1, k86=1, k34=1, k43=1, k56=1, k65=1, 
    k12=1, k21=1, k78=1, k87=1, k71=1, k17=1, 
    k53=1, k35=1, k64=1, k46=1, k82=1, k28=1,
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


def get_K_vals(k_AA_anti, k_AA_sym):
    H_on = 1e10
    H_off = 1e3
    D_on = 1e7
    D_off = 1
    H_in = 10 ** (-6.5)
    H_out = 10 ** (-7.5)
    c_D = 25e-9

    K = get_K(
        k31=H_on * H_out,
        k13=H_off,
        k57=H_on * H_in,
        k75=H_off,
        k42=H_on * H_out,
        k24=H_off,
        k68=H_on * H_in,
        k86=H_off,
        k34=D_on * c_D,
        k43=D_off,
        k56=D_on * c_D,
        k65=D_off,
        k12=D_on * c_D,
        k21=D_off,
        k78=D_on * c_D,
        k87=D_off,
        k71=k_AA_anti,
        k17=k_AA_anti,
        k53=k_AA_sym,
        k35=k_AA_sym,
        k64=k_AA_anti,
        k46=k_AA_anti,
        k82=k_AA_sym,
        k28=k_AA_sym,
    )
    return K

def get_K_and_G():
    # generate array using defaults (1)
    K = get_K()
    G = nx.MultiDiGraph()
    graph_utils.generate_edges(G, K)
    return K, G


def get_rate_names(fig_key):
    if fig_key == "7A":
        return ["k_AA_anti", "k_AA_sym"]
    elif fig_key == "7B":
        return ["k_EH_H", "k_EHD_H", "k_ED_D", "k_EHD_D", "k_AA"]


def get_sub_dict(fig_key):
    (k31,k13,k57,k75,k42,k24,k68,k86,k34,k43,k56,k65,k12,k21,k78,k87,k71,k17,k53,k35,k64,k46,k82,k28) = symbols(
        "k31 k13 k57 k75 k42 k24 k68 k86 k34 k43 k56 k65 k12 k21 k78 k87 k71 k17 k53 k35 k64 k46 k82 k28")

    H_on = 1e10
    H_off = 1e3
    D_on = 1e7
    D_off = 1
    H_in = 10 ** (-6.5)
    H_out = 10 ** (-7.5)
    c_D = 25e-9

    if fig_key == "7A":
        (k_AA_anti, k_AA_sym) = symbols("k_AA_anti k_AA_sym")
        sub_dict = {
            k31: H_on * H_out,
            k13: H_off,
            k57: H_on * H_in,
            k75: H_off,
            k42: H_on * H_out,
            k24: H_off,
            k68: H_on * H_in,
            k86: H_off,
            k34: D_on * c_D,
            k43: D_off,
            k56: D_on * c_D,
            k65: D_off,
            k12: D_on * c_D,
            k21: D_off,
            k78: D_on * c_D,
            k87: D_off,
            k71: k_AA_anti,
            k17: k_AA_anti,
            k53: k_AA_sym,
            k35: k_AA_sym,
            k64: k_AA_anti,
            k46: k_AA_anti,
            k82: k_AA_sym,
            k28: k_AA_sym,
        }
    elif fig_key == "7B":
        (k_EH_H, k_EHD_H, k_ED_D, k_EHD_D, k_AA) = symbols(
            "k_EH_H k_EHD_H k_ED_D k_EHD_D k_AA")
        sub_dict = {
            k31: H_on * H_out,
            k13: k_EH_H,
            k57: H_on * H_in,
            k75: k_EH_H,
            k42: H_on * H_out,
            k24: k_EHD_H,
            k68: H_on * H_in,
            k86: k_EHD_H,
            k34: D_on * c_D,
            k43: k_ED_D,
            k56: D_on * c_D,
            k65: k_ED_D,
            k12: D_on * c_D,
            k21: k_EHD_D,
            k78: D_on * c_D,
            k87: k_EHD_D,
            k71: k_AA,
            k17: k_AA,
            k53: k_AA,
            k35: k_AA,
            k64: k_AA,
            k46: k_AA,
            k82: k_AA,
            k28: k_AA,
        }
    else:
        raise ValueError("Invalid option, please choose 7A or 7B")


    return sub_dict


def generate_probability_expressions(G):

    # all figures share the probability expressions
    expression_path = "probability_expressions"
    expr_filepath = os.path.abspath(f"{expression_path}/p_1.txt")

    if os.path.isfile(expr_filepath):
        # if the probability expressions already exist, load them
        # in and append them to a list
        print("--> Probability expressions found, collecting now...")
        sys.stdout.flush()

        pi_expressions = []
        for i in range(G.number_of_nodes()):
            _fpath = os.path.abspath(f"{expression_path}/p_{i+1}.txt")
            with open(_fpath, "r") as f:
                _expr = f.readline()
                pi_expressions.append(_expr)
    else:
        # if the expressions do not exist, create them and save them

        # create the transition flux expressions using generic figure
        print("--> Probability expressions not found. Generating expressions...")
        sys.stdout.flush()

        pi_expressions = calculations.calc_state_probs(G, key="name", output_strings=True)

        print(f"--> Expression generation complete. Storing expressions in ~/{expression_path}")
        for i, _expr in enumerate(pi_expressions):
            _fpath = os.path.abspath(f"{expression_path}/p_{i+1}.txt")
            with open(_fpath, "w") as f:
                f.write(str(_expr))

    return pi_expressions


def generate_sympy_funcs(expressions, fig_key):

    func_path = "sympy_funcs"

    # get the dictionary of variable substitutions
    sub_dict = get_sub_dict(fig_key=fig_key)

    # substitute variables and simplify the expressions
    sympy_funcs = []
    for i, func in enumerate(expressions):
        _fpath = os.path.abspath(f"{func_path}/p_{i+1}_{fig_key}.txt")
        if not os.path.isfile(_fpath):
            print(
                f"--> No SymPy function found at location {_fpath} \n"
                f"--> Generating new SymPy function..."
            )
            func = parse_expr(func)
            simplified_func = simplify(func.subs(sub_dict))
            dill.dump(simplified_func, open(_fpath, "wb"))
        else:
            print(f"--> Loading SymPy function from location {_fpath}")
            simplified_func = dill.load(open(_fpath, "rb"))
        sys.stdout.flush()
        sympy_funcs.append(simplified_func)

    return sympy_funcs


def generate_lambda_funcs(sympy_funcs, fig_key):
    # collect the rate names for the order parameters
    rate_names = get_rate_names(fig_key=fig_key)

    lambda_funcs = []
    for _func in sympy_funcs:
        _func_lambda = lambdify(rate_names, _func, "numpy")
        lambda_funcs.append(_func_lambda)

    return lambda_funcs
