import sys
# have to raise the recursion limit or SymPy raises a RecursionError when
# attempting to parse the EmrE 8-state model expression for Sigma
sys.setrecursionlimit(5000)
import numpy as np
import networkx as nx
from sympy import symbols
from sympy.parsing.sympy_parser import parse_expr

from kda import graph_utils, diagrams, calculations as calcs


# define connectivity matrix
K = np.array(
    [
        [0, 1, 1, 0, 0, 0, 1, 0],
        [1, 0, 0, 1, 0, 0, 0, 1],
        [1, 0, 0, 1, 1, 0, 0, 0],
        [0, 1, 1, 0, 0, 1, 0, 0],
        [0, 0, 1, 0, 0, 1, 1, 0],
        [0, 0, 0, 1, 1, 0, 0, 1],
        [1, 0, 0, 0, 1, 0, 0, 1],
        [0, 1, 0, 0, 0, 1, 1, 0],
    ]
)
# create the 8-state kinetic diagram
G = nx.MultiDiGraph()
# populate edge data using KDA utility
graph_utils.generate_edges(G, K)
# create symbols for current variables
(k31, k13, k57, k75, k42, k24, k68, k86, k34, k43, k56, k65, 
	k12, k21, k78, k87, k71, k17, k53, k35, k64, k46, k82, k28) = symbols(
	"k31 k13 k57 k75 k42 k24 k68 k86 k34 k43 k56 k65 k12 k21 k78 k87 k71 k17 k53 k35 k64 k46 k82 k28")
# create variables for substitution
(H_on, H_off, D_on, D_off, H_in, H_out, D_in, D_out, k_AA_anti, k_AA_sym) = symbols(
    "H_on H_off D_on D_off H_in H_out D_in D_out k_AA_anti k_AA_sym")
# define variable substitution mapping
model = {
    k31: H_on * H_out, k13: H_off, 
    k57: H_on * H_in, k75: H_off, 
    k42: H_on * H_out, k24: H_off, 
    k68: H_on * H_in, k86: H_off, 
    k34: D_on * D_out, k43: D_off, 
    k56: D_on * D_in, k65: D_off, 
    k12: D_on * D_out, k21: D_off, 
    k78: D_on * D_in, k87: D_off, 
    k71: k_AA_anti, k17: k_AA_anti, 
    k53: k_AA_sym, k35: k_AA_sym, 
    k64: k_AA_anti, k46: k_AA_anti, 
    k82: k_AA_sym, k28: k_AA_sym,
}

# collect the directional edges beforehand
dir_edges = diagrams.generate_directional_diagrams(G, return_edges=True)
# create the normalization factor for net cycle fluxes
sigma = calcs.calc_sigma(G, dir_edges, key="name", output_strings=True)
sigma = parse_expr(sigma)
sigma = sigma.subs(model).simplify()

# define cycle info
EmrE_cycles = {
	1: {"cycle": [0, 1, 3, 2, 4, 5, 7, 6], "order": [6, 0]}, 
	2: {"cycle": [0, 6, 7, 5, 4, 2], "order": [6, 0]}, 
	3: {"cycle": [0, 6, 7, 5, 3, 2], "order": [6, 0]}, 
	4: {"cycle": [0, 6, 7, 5, 3, 1], "order": [6, 0]}, 
	6: {"cycle": [0, 6, 7, 1, 3, 2], "order": [6, 0]}, 
	7: {"cycle": [0, 6, 7, 1], "order": [6, 0]}, 
	8: {"cycle": [0, 2, 3, 1, 7, 5, 4, 6], "order": [6, 0]}, 
	9: {"cycle": [0, 6, 4, 5, 7, 1], "order": [6, 0]}, 
	10: {"cycle": [0, 6, 4, 5, 3, 2], "order": [6, 0]}, 
	11: {"cycle": [0, 6, 4, 5, 3, 1], "order": [6, 0]}, 
	13: {"cycle": [0, 6, 4, 2, 3, 1], "order": [6, 0]}, 
	14: {"cycle": [0, 6, 4, 2], "order": [6, 0]}, 
	15: {"cycle": [0, 1, 3, 5, 7, 6, 4, 2], "order": [0, 2]}, 
	16: {"cycle": [0, 2, 4, 6, 7, 1], "order": [0, 2]}, 
	17: {"cycle": [0, 2, 4, 5, 7, 1], "order": [0, 2]}, 
	18: {"cycle": [0, 2, 4, 5, 3, 1], "order": [0, 2]}, 
	19: {"cycle": [0, 2, 3, 5, 7, 1], "order": [0, 2]}, 
	20: {"cycle": [0, 1, 7, 6, 4, 5, 3, 2], "order": [0, 2]}, 
	22: {"cycle": [1, 7, 6, 4, 5, 3], "order": [1, 3]}, 
	23: {"cycle": [1, 7, 6, 4, 2, 3], "order": [2, 4]}, 
	24: {"cycle": [1, 7, 5, 4, 2, 3], "order": [2, 4]}, 
	25: {"cycle": [1, 7, 5, 3], "order": [1, 3]}, 
	26: {"cycle": [2, 3, 5, 7, 6, 4], "order": [2, 4]}, 
	27: {"cycle": [2, 3, 5, 4], "order": [2, 4]},
	}

# generate net cycle flux expression numerators
for idx, cycle_info in EmrE_cycles.items():
	pi_diff = calcs.calc_pi_difference(G, cycle=cycle_info["cycle"], 
		order=cycle_info["order"], key="name", output_strings=True)
	flux_diags = diagrams.generate_flux_diagrams(G, cycle=cycle_info["cycle"])
	sigma_K = calcs.calc_sigma_K(G, cycle=cycle_info["cycle"],
		flux_diags=flux_diags, key="name", output_strings=True)
	if sigma_K == 1:
		func = parse_expr(pi_diff)
	else:
		func = parse_expr(pi_diff) * parse_expr(sigma_K)
	# perform variable substitutions
	func = func.subs(model).simplify()
	EmrE_cycles[idx]["func"] = func


# define contributing cycles for both ligands
H_cycles = [1, 2, 3, 4, 10, 11, 13, 14, 16, 17, 19, 20, 22, 23, 24, 25]
D_cycles = [3, 4, 6, 7, 8, 9, 10, 11, 15, 16, 17, 18, 23, 24, 26, 27]

J_H_cycle = 0
for idx in H_cycles:
	J_H_cycle += EmrE_cycles[idx]["func"]

J_D_cycle = 0
for idx in D_cycles:
	# negative contributors
	if idx in (3, 4, 6, 7, 8, 9, 10, 11):
		J_D_cycle -= EmrE_cycles[idx]["func"]
	else:
		J_D_cycle += EmrE_cycles[idx]["func"]

# normalize the operational fluxes (cycles)
J_H_cycle = (J_H_cycle / sigma).simplify()
J_D_cycle = (J_D_cycle / sigma).simplify()

# generate state probability expressions
prob_strs = calcs.calc_state_probs(G, key='name', output_strings=True)
p1, p2, p3, p4, p5, p6, p7, p8 = prob_strs

# calculate operational fluxes using transition fluxes
J_13 = p1.subs(model) * model[k13] - p3.subs(model) * model[k31]
J_24 = p2.subs(model) * model[k24] - p4.subs(model) * model[k42]
J_21 = p2.subs(model) * model[k21] - p1.subs(model) * model[k12]
J_43 = p4.subs(model) * model[k43] - p3.subs(model) * model[k34]
J_H_trans = (J_13 + J_24).simplify()
J_D_trans = (J_21 + J_43).simplify()

# compare operational flux expressions
# created using cycles and transitions
print((J_H_cycle - J_H_trans).simplify() == 0)
print((J_D_cycle - J_D_trans).simplify() == 0)

# =================
# === Output ======
# =================
# (kda-env)
# nikol@Lorentz MINGW64 ~/OneDrive/projects/kinetic_diagram_analysis/kda-examples/kda_analysis/emre (add_op_flux_comparison_scripts)
# $ python emre_operational_flux_expression_comparison.py
# Cycle [0, 1, 3, 2, 4, 5, 7, 6] contains all nodes in G. No flux diagrams generated.
# No flux diagrams detected for cycle [0, 1, 3, 2, 4, 5, 7, 6]. Sigma K value is 1.
# Cycle [0, 2, 3, 1, 7, 5, 4, 6] contains all nodes in G. No flux diagrams generated.
# No flux diagrams detected for cycle [0, 2, 3, 1, 7, 5, 4, 6]. Sigma K value is 1.
# Cycle [0, 1, 3, 5, 7, 6, 4, 2] contains all nodes in G. No flux diagrams generated.
# No flux diagrams detected for cycle [0, 1, 3, 5, 7, 6, 4, 2]. Sigma K value is 1.
# Cycle [0, 1, 7, 6, 4, 5, 3, 2] contains all nodes in G. No flux diagrams generated.
# No flux diagrams detected for cycle [0, 1, 7, 6, 4, 5, 3, 2]. Sigma K value is 1.
# True
# True