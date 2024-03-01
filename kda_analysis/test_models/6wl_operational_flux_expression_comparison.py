import numpy as np
import networkx as nx
from sympy import symbols

from kda import graph_utils, calculations as calcs


# define connectivity matrix
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
# create the 6-state kinetic diagram
G_leak = nx.MultiDiGraph()
# populate edge data using KDA utility
graph_utils.generate_edges(G_leak, K)

# create symbols for current variables
(k12, k21, k23, k32, k34, k43, k45, 
    k54, k56, k65, k61, k16, k14, k41) = symbols(
    "k12 k21 k23 k32 k34 k43 k45 k54 k56 k65 k61 k16 k14, k41")

# create variables for substitution
(H_on, H_off, Na_on, Na_off, k_conf, 
k_leak,  H_in, H_out, c_Na) = symbols(
    "H_on, H_off, Na_on, Na_off, k_conf, k_leak,  H_in, H_out, c_Na")

# define variable substitution mapping
model = {
    k12:H_on*H_out, k21:H_off,
    k23:k_conf, k32:k_conf,
    k34:H_off, k43:H_on*H_in,
    k45:Na_on*c_Na, k54:Na_off,
    k56:k_conf, k65:k_conf,
    k61:Na_off, k16:Na_on*c_Na,
    k14:k_leak, k41:k_leak,
}

# collect the cycles using KDA
cycles = graph_utils.find_all_unique_cycles(G_leak)
# label cycles according to figure
cycle_labels = ["a", "c", "b"]
# set the cycle directions (CW)
cycle_orders = [[0, 1], [0, 1], [5, 0]]

# generate net cycle flux expressions
G_leak_cycles = {}
for label, _cycle, _order in zip(cycle_labels, cycles, cycle_orders):
    func = calcs.calc_net_cycle_flux(
        G_leak, cycle=_cycle, order=_order, 
        key='name', output_strings=True,
        )
    # perform variable substitutions
    func = func.subs(model).simplify()
    G_leak_cycles[label] = {
        "cycle": _cycle, 
        "order": _order, 
        "func": func,
        }

# assign net cycle flux expressions
J_a = G_leak_cycles["a"]["func"]
J_b = G_leak_cycles["b"]["func"]
J_c = G_leak_cycles["c"]["func"]

# generate state probability expressions
prob_strs = calcs.calc_state_probs(
    G_leak, key='name', output_strings=True)
p1, p2, p3, p4, p5, p6 = prob_strs

# calculate operational fluxes using cycles
J_H_cycle = (J_a + J_c).simplify()
J_Na_cycle = (J_b + J_c).simplify()

# calculate operational fluxes using transition fluxes
# J_H = p1 * k12 - p2 * k21 = J12
J_H_trans = (p1.subs(model) * model[k12] 
    - p2.subs(model) * model[k21]).simplify()
# J_Na = p6 * k61 - p1 * k16 = J61
J_Na_trans = (p6.subs(model) * model[k61]
    - p1.subs(model) * model[k16]).simplify()

# compare operational flux expressions
# created using cycles and transitions
print((J_H_cycle - J_H_trans).simplify() == 0)
print((J_Na_cycle - J_Na_trans).simplify() == 0)

# =================
# === Output ======
# =================
# (kda-env)
# nikol@Lorentz MINGW64 ~/OneDrive/projects/kinetic_diagram_analysis/kda-examples/kda_analysis/test_models (add_op_flux_comparison_scripts)
# $ python 6wl_operational_flux_expression_comparison.py
# Cycle [0, 1, 2, 3, 4, 5] contains all nodes in G. No flux diagrams generated.
# No flux diagrams detected for cycle [0, 1, 2, 3, 4, 5]. Sigma K value is 1.
# True
# True