Kinetic Diagram Analysis Examples
====================================
[//]: # (Badges)
![CI](https://github.com/Becksteinlab/kda-examples/actions/workflows/main_ci.yml/badge.svg)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.10795714.svg)](https://doi.org/10.5281/zenodo.10795714)

Example code for building models using [Kinetic Diagram Analysis](https://github.com/Becksteinlab/kda).

## Examples

<img src="kda_examples/test_model_3_state/diagrams/input.png"  width=300 alt="test model 3-state"> <img src="kda_examples/test_model_4_state/diagrams/input.png" width=300, alt="test model 4-state">
<img src="kda_examples/test_model_4_state_leakage/diagrams/input.png" width=300, alt="test model 4-state leakage"> <img src="kda_examples/test_model_5_state_leakage/diagrams/input.png" width=300, alt="test model 5-state leakage">
<img src="kda_examples/test_model_6_state/diagrams/input.png" width=300, alt="test model 6-state"> <img src="kda_examples/test_model_6_state_leakage/diagrams/input.png" width=300, alt="test model 6-state leakage">
<img src="kda_examples/test_model_8_state_leakage/diagrams/input.png" width=300, alt="test model 8-state leakage">

## Adding Examples

See [PR #2](https://github.com/Becksteinlab/kda-examples/pull/2) for example on the minimum requirements for adding new examples to the repo.

## Testing

This package has a testing module that checks if the example scripts execute
without error. To run this test, simply install [KDA](https://github.com/Becksteinlab/kda), install
this package (instructions below) and run `pytest`.

## Installation
### Development version from source

To install the latest development version from source, run
```bash
git clone git@github.com:Becksteinlab/kda-examples.git
cd kda-examples
python setup.py install
```

## Copyright

Copyright (c) 2022, Nikolaus Awtrey
