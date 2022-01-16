import os
import numpy as np
from kda_examples.generate_diagrams import generate_diagrams


def main():
    # define the input/connectivity matrix
    K = np.array(
        [
            [0, 1, 0, 1],
            [1, 0, 1, 0],
            [0, 1, 0, 1],
            [1, 0, 1, 0],
        ]
    )
    # specify the positions of nodes in NetworkX fashion
    node_positions = {0: [1, 1], 1: [-1, 1], 2: [-1, -1], 3: [1, -1]}
    # get path to /diagrams
    save_path = os.path.abspath("diagrams")
    # generate partial, directional, and flux diagrams
    # and save in the /diagrams directory
    generate_diagrams(
        input_mat=K,
        node_positions=node_positions,
        save_path=save_path,
        # NOTE: there are no flux diagrams for the 4-state model
        gen_flux_diagrams=False
    )

if __name__ == "__main__":
    main()
