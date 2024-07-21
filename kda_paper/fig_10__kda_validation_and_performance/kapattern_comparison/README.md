KDA vs. KAPattern Expression Comparison
=======================================

The included script, `comparison.py`, was written with the purpose of comparing
`KDA` expressions to the expressions output by the `KAPattern` software using 
`sympy`. The script uses previously pickled `NetworkX.MultiDiGraph` objects 
created by the `KDA` [verification script](https://github.com/Becksteinlab/kda/tree/master/kda/scripts). The specific graphs included 
here are from the data set listed in the [`kda-examples`](https://github.com/Becksteinlab/kda-examples/tree/master/kda_paper/fig_10__kda_validation_and_performance) repository.

For reproducibility purposes, our pickled graph objects, input files, and
output files are included here. The terminal output from our comparison run
is included in `comparison_output.txt`.

## Requirements

- `KDA v0.3.0` (see [installation instructions](https://github.com/becksteinlab/kda?tab=readme-ov-file#installation))
- `KAPattern v1.1` (download [here](https://vpr.sites.uofmhosting.net/software/ka-pattern))


## Comparing Expressions

To compare `KDA` and `KAPattern` expressions the following steps
must be completed:

1. Initial `comparison.py` run: generates the `KAPattern` input files
   (e.g. `input_graph_3_1.txt`), plots and saves the graphs in the current
   working directory as `.pdf` files
2. Using `KAPattern` and the generated input files, individually run and
   save the output files in the same format as the input files 
   (e.g. `output_graph_3_1.txt`)
3. Final `comparison.py` run: outputs the comparison results in the terminal


### Running `comparison.py`

After activating your `KDA` environment, run the script in the same
directory as the `.pk` files:

```bash
1  cd ~/path/to/directory
2  python comparison.py
```

The script is setup to run the following cases from the KDA Validation Data:
- `graph_3_1.pk`
- `graph_4_1.pk`
- `graph_5_1.pk`
- `graph_6_9.pk`
- `graph_7_17.pk`


### `comparison.py` Outputs

The script will plot and save the graphs in the current working directory
as `.pdf` files and print out the results of the comparison for each model:

```bash
Loaded graph_3_1.pk -- MultiDiGraph with 3 nodes and 6 edges

KDA vs KAPattern State Probability Expression Comparison
========================================
File: graph_3_1.pk
Graph: MultiDiGraph with 3 nodes and 6 edges
----------------------------------------
>> State 1 expressions agree: True
>> State 2 expressions agree: True
>> State 3 expressions agree: True
========================================
```

If expressions do not agree, the resultant simplified expression will be
shown (to highlight the disagreement) and the full expressions will be
displayed for comparison:

```bash
Loaded graph_6_9.pk -- MultiDiGraph with 6 nodes and 18 edges

KDA vs KAPattern State Probability Expression Comparison
========================================
File: graph_6_9.pk
Graph: MultiDiGraph with 6 nodes and 18 edges
----------------------------------------
>> State 1 expressions agree: False
>>>> Expressions do not agree! The resultant expression: 
	-k25*k34*k41*k45*k54*k63 + k25*k34*k41*k54*k63
>>>> KDA Expression: 
	k21*k31*k41*k52*k63 + k21*k31*k41*k52*k65 + k21*k31*k41*k54*k63 
	+ k21*k31*k41*k54*k65 + k21*k31*k41*k56*k63 + k21*k31*k42*k52*k63 
	+ k21*k31*k42*k52*k65 + k21*k31*k42*k54*k63 + k21*k31*k42*k54*k65 
	+ k21*k31*k42*k56*k63 + k21*k31*k43*k52*k63 + k21*k31*k43*k52*k65 
	+ k21*k31*k43*k54*k63 + k21*k31*k43*k54*k65 + k21*k31*k43*k56*k63 
	+ k21*k31*k45*k52*k63 + k21*k31*k45*k52*k65 + k21*k31*k45*k56*k63 
	+ k21*k34*k41*k52*k63 + k21*k34*k41*k52*k65 + k21*k34*k41*k54*k63 
	+ k21*k34*k41*k54*k65 + k21*k34*k41*k56*k63 + k21*k34*k42*k52*k63 
	+ k21*k34*k42*k52*k65 + k21*k34*k42*k54*k63 + k21*k34*k42*k54*k65 
	+ k21*k34*k42*k56*k63 + k21*k34*k45*k52*k63 + k21*k34*k45*k52*k65 
	+ k21*k36*k41*k52*k65 + k21*k36*k41*k54*k65 + k21*k36*k42*k52*k65 
	+ k21*k36*k42*k54*k65 + k21*k36*k43*k52*k65 + k21*k36*k45*k52*k65 
	+ k24*k31*k41*k52*k63 + k24*k31*k41*k52*k65 + k24*k31*k41*k54*k63 
	+ k24*k31*k41*k54*k65 + k24*k31*k41*k56*k63 + k24*k31*k43*k52*k63 
	+ k24*k31*k43*k52*k65 + k24*k31*k43*k54*k63 + k24*k31*k43*k54*k65 
	+ k24*k31*k43*k56*k63 + k24*k31*k45*k56*k63 + k24*k34*k41*k52*k63 
	+ k24*k34*k41*k52*k65 + k24*k34*k41*k54*k63 + k24*k34*k41*k54*k65 
	+ k24*k34*k41*k56*k63 + k24*k36*k41*k52*k65 + k24*k36*k41*k54*k65 
	+ k25*k31*k41*k54*k63 + k25*k31*k41*k54*k65 + k25*k31*k41*k56*k63 
	+ k25*k31*k42*k56*k63 + k25*k31*k43*k54*k63 + k25*k31*k43*k54*k65 
	+ k25*k31*k43*k56*k63 + k25*k31*k45*k56*k63 + k25*k34*k41*k54*k63 
	+ k25*k34*k41*k54*k65 + k25*k34*k41*k56*k63 + k25*k36*k41*k54*k65
>>>> KAPattern Expression: 
	k21*k31*k41*k52*k63 + k21*k31*k41*k52*k65 + k21*k31*k41*k54*k63 
	+ k21*k31*k41*k54*k65 + k21*k31*k41*k56*k63 + k21*k31*k42*k52*k63 
	+ k21*k31*k42*k52*k65 + k21*k31*k42*k54*k63 + k21*k31*k42*k54*k65 
	+ k21*k31*k42*k56*k63 + k21*k31*k43*k52*k63 + k21*k31*k43*k52*k65 
	+ k21*k31*k43*k54*k63 + k21*k31*k43*k54*k65 + k21*k31*k43*k56*k63 
	+ k21*k31*k45*k52*k63 + k21*k31*k45*k52*k65 + k21*k31*k45*k56*k63 
	+ k21*k34*k41*k52*k63 + k21*k34*k41*k52*k65 + k21*k34*k41*k54*k63 
	+ k21*k34*k41*k54*k65 + k21*k34*k41*k56*k63 + k21*k34*k42*k52*k63 
	+ k21*k34*k42*k52*k65 + k21*k34*k42*k54*k63 + k21*k34*k42*k54*k65 
	+ k21*k34*k42*k56*k63 + k21*k34*k45*k52*k63 + k21*k34*k45*k52*k65 
	+ k21*k36*k41*k52*k65 + k21*k36*k41*k54*k65 + k21*k36*k42*k52*k65 
	+ k21*k36*k42*k54*k65 + k21*k36*k43*k52*k65 + k21*k36*k45*k52*k65 
	+ k24*k31*k41*k52*k63 + k24*k31*k41*k52*k65 + k24*k31*k41*k54*k63 
	+ k24*k31*k41*k54*k65 + k24*k31*k41*k56*k63 + k24*k31*k43*k52*k63 
	+ k24*k31*k43*k52*k65 + k24*k31*k43*k54*k63 + k24*k31*k43*k54*k65 
	+ k24*k31*k43*k56*k63 + k24*k31*k45*k56*k63 + k24*k34*k41*k52*k63 
	+ k24*k34*k41*k52*k65 + k24*k34*k41*k54*k63 + k24*k34*k41*k54*k65 
	+ k24*k34*k41*k56*k63 + k24*k36*k41*k52*k65 + k24*k36*k41*k54*k65 
	+ k25*k31*k41*k54*k63 + k25*k31*k41*k54*k65 + k25*k31*k41*k56*k63 
	+ k25*k31*k42*k56*k63 + k25*k31*k43*k54*k63 + k25*k31*k43*k54*k65 
	+ k25*k31*k43*k56*k63 + k25*k31*k45*k56*k63 + k25*k34*k41*k45*k54*k63 
	+ k25*k34*k41*k54*k65 + k25*k34*k41*k56*k63 + k25*k36*k41*k54*k65
```