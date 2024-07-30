"""
File for testing the evaluation time of KDA expressions
using the timeit module for some of the randomly
generated graphs from the KDA Verification data.
"""

import os
import timeit
import pickle
from sympy import lambdify

from kda import calculations


def load_pickled_graph(path):
	"""
	Loads a pickled `NetworkX.MultiDiGraph` object
	"""
	with open(path, 'rb') as f:
		G = pickle.load(f)
	return G


def test(func, vals):
	return func(**vals)


if __name__ == "__main__":
	graph_paths = [
		os.path.abspath("./graph_3_1.pk"),
		os.path.abspath("./graph_4_1.pk"),
		os.path.abspath("./graph_5_1.pk"),
		os.path.abspath("./graph_6_9.pk"),
		os.path.abspath("./graph_7_17.pk"),
	]

	graphs = []
	times = []
	for path in graph_paths:

		G = load_pickled_graph(path=path)
		p1 = calculations.calc_state_probs(G=G, key="name", output_strings=True)[0]

		var_names = list(p1.free_symbols)
		var_names = [str(v) for v in var_names]
		vals = {name:i+1 for i, name in enumerate(var_names)}
		p1_lambda = lambdify(var_names, p1, "numpy")

		p1_val = p1_lambda(**vals)
		print(f"p1 for {os.path.basename(path)}: {p1_val}")

		iterations = 1000
		t_avg = timeit.timeit("test(func=p1_lambda, vals=vals)", number=iterations, globals=globals())/iterations
		graphs.append(os.path.basename(path))
		times.append(t_avg)

	print(f"Average runtimes for {iterations} iterations")
	print("-" * 40)
	for g, t in zip(graphs, times):
		print(f"{g}:  {t:.1e}")


