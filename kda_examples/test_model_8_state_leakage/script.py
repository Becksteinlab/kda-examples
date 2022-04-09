import os
import numpy as np
from kda_examples.generate_diagrams import generate_diagrams


def main():
    # define the input/connectivity matrix
    K = np.array(
        [
            [0, 1, 1, 0, 0, 0, 0, 0],
            [1, 0, 0, 1, 0, 0, 0, 0],
            [1, 0, 0, 1, 1, 0, 0, 0],
            [0, 1, 1, 0, 0, 1, 0, 0],
            [0, 0, 1, 0, 0, 1, 1, 0],
            [0, 0, 0, 1, 1, 0, 0, 1],
            [0, 0, 0, 0, 1, 0, 0, 1],
            [0, 0, 0, 0, 0, 1, 1, 0],
        ],
    )
    # specify the positions of nodes in NetworkX fashion
    node_positions = {
        0 : [-0.5, 1.5],
        1 : [0.5, 1.5],
        2 : [-0.5, 0.5],
        3 : [0.5, 0.5],
        4 : [-0.5, -0.5],
        5 : [0.5, -0.5],
        6 : [-0.5, -1.5],
        7 : [0.5, -1.5],
    }
    # get path to /diagrams
    save_path = os.path.abspath("diagrams")
    # generate partial, directional, and flux diagrams
    # and save in the /diagrams directory
    generate_diagrams(
        input_mat=K,
        node_positions=node_positions,
        save_path=save_path,
        gen_flux_diagrams=True,
        set_panel_rows=False,
    )


if __name__ == "__main__":
    main()
