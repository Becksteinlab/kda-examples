import sys
# raise the recursion limit so SymPy doesn't raise a
# RecursionError when parsing larger expressions
sys.setrecursionlimit(5000)
import os

import pickle
from sympy import fraction
from sympy.parsing.sympy_parser import parse_expr

from kda import expressions, calculations, plotting


def load_pickled_graph(path):
	"""
	Loads a pickled `NetworkX.MultiDiGraph` object
	"""
	with open(path, 'rb') as f:
		G = pickle.load(f)
	print(f"Loaded {os.path.basename(path)} -- {G}")
	return G


def get_KAPattern_input_filename(path):
	"""
	Creates a new file path of the form `~/input_graph_X_X.txt`
	from an input absolute file path.
	"""
	# use the absolute path to create a new file name
	# in the format `./input_graph_3_1.txt`
	fname = f"input_{os.path.splitext(os.path.basename(path))[0]}.txt"
	return os.path.abspath(fname)


def get_KAPattern_output_filename(path):
	"""
	Creates a new file path of the form `~/output_graph_X_X.txt`
	from an input absolute file path.
	"""
	fname = f"output_{os.path.splitext(os.path.basename(path))[0]}.txt"
	return os.path.abspath(fname)


def write_KAPattern_input_file_from_graph(G, save_path):
	"""
	Takes a `NetworkX.MultiDiGraph` and stores its
	information in KAPattern format with generic
	edge labels (i.e. `kij`).
	"""
	with open(save_path, 'w') as f:
		print(f"Writing graph data in {save_path}")
		for i, j, attrs in G.edges(data=True):
			f.write(f"{i+1} {j+1} {attrs['name']}\n")


def read_KAPattern_output_file(output_path):
	"""
	Opens a `KAPattern` output file and retrieves the
	state multiplicity expression strings
	"""
	exprs = []
	with open(output_path, "r") as f:
		for line in f:
			# E(1) lines always start with a space
			if (line.startswith(" E(1)") or line.startswith("E(")):
				# for each line, grab the expression after the equal sign
				# and store it. Since each expression is in increasing
				# consecutive order they should be stored in the correct
				# order in the list
				expr = line.split(" = ")[1].split(";")[0]
				exprs.append(expr)
	return exprs



def main():
	graph_paths = [
		os.path.abspath("./graph_3_1.pk"),
		os.path.abspath("./graph_4_1.pk"),
		os.path.abspath("./graph_5_1.pk"),
		os.path.abspath("./graph_6_9.pk"),
		os.path.abspath("./graph_7_17.pk"),
	]
	for path in graph_paths:
		G = load_pickled_graph(path=path)

		input_fpath = get_KAPattern_input_filename(path=path)
		output_fpath = get_KAPattern_output_filename(path=path)

		if not os.path.isfile(output_fpath):
			# plot the graph and save it
			fig = plotting.draw_diagrams(G)
			fig_fname = os.path.splitext(path)[0] + ".pdf"
			fig.savefig(fig_fname, dpi=300)
			print("No output files detected. Generating input file...")
			write_KAPattern_input_file_from_graph(G=G, save_path=input_fpath)
			msg = (
					f"Please re-run this script after 1.) running input files through `KAPattern.exe`"
					f" and 2.) saving the output files in the format `output_graph_3_1.txt`."
				)
			print(msg)
		else:
			mult_exprs = read_KAPattern_output_file(output_path=output_fpath)
			# convert strings to sympy expressions
			KAP_prob_exprs = [parse_expr(e) for e in mult_exprs]

			KDA_prob_exprs = calculations.calc_state_probs(G=G, key="name", output_strings=True)
			# retrieve the numerator from each normalized state probability expression
			KDA_prob_exprs = [fraction(e)[0] for e in KDA_prob_exprs]

			print(f"\nKDA vs KAPattern State Probability Expression Comparison")
			print("=" * 40)
			print(f"File: {os.path.basename(path)}")
			print(f"Graph: {G}")
			print("-" * 40)
			for i, (KDA_expr, KAP_expr) in enumerate(zip(KDA_prob_exprs, KAP_prob_exprs)):
				# see if the expressions will fully cancel
				comparison_expr = (KDA_expr - KAP_expr).cancel()
				truthVal = comparison_expr == 0
				print(f">> State {i+1} expressions agree: {truthVal}")
				if not truthVal:
					print(f">>>> Expressions do not agree! The resultant expression: \n{comparison_expr}")
					print(f">>>> KDA Expression: \n{KDA_expr}")
					print(f">>>> KAPattern Expression: \n{KAP_expr}")
			print("=" * 40 + "\n\n")
			sys.stdout.flush()


if __name__ == "__main__":
	main()
