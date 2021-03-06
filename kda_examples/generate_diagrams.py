import networkx as nx
from kda import graph_utils, diagrams, plotting


def _flatten_flux_diagrams(diagrams):
    # temporary work around to flatten output of
    # `kda.diagrams.generate_all_flux_diagrams()`
    flux_diagrams = []
    for diag_list in diagrams:
        if diag_list:
            for diag in diag_list:
                flux_diagrams.append(diag)
    return flux_diagrams


def generate_diagrams(input_mat, node_positions, save_path, gen_flux_diagrams=True, set_panel_rows=True):
    """
    Constructs a kinetic diagram based on the input matrix, generates the
    partial, directional, and flux diagrams for that diagram, and stores
    the outputs in the `/diagrams` directory.

    Parameters
    ----------
    input_mat: array
        'NxN' array where 'N' is the number of nodes.

    node_positions: dict
        Dictionary where keys are the indexed states (0, 1, 2, ..., N)
        and the values are the x, y coordinates for each node.

    save_path: str
        String of save path for figure. If a path is specified the figure
        will be saved at the specified location.

    gen_flux_diagrams: bool
        Binary used to conveniently switch on/off the flux diagram generation
        since some models do not have any flux diagrams.

    set_panel_rows: bool
        Binary used to determine whether partial and directional diagram
        panels have 1 row and `n_states` rows, respectively (`True`), or
        use the `draw_diagrams()` built-in method for determining panel
        dimensions (`False`). Default is `True`.

    """
    if set_panel_rows:
        # set partial diagram panel rows to 1
        partial_rows = 1
        # set directional diagram panel rows to number of states
        directional_rows = input_mat.shape[0]
    else:
        # pass `None` so the panel dimensions are set by
        # `plotting.draw_diagrams()`. This will generate
        # a panel that is roughly square
        partial_rows = None
        directional_rows = None

    # initialize an empty graph object
    G = nx.MultiDiGraph()
    # populate the edge data
    graph_utils.generate_edges(G, input_mat)

    #==============================
    #=== Generate the diagrams ====
    #==============================
    # generate the partial and directional diagrams for G
    partial_diagrams = diagrams.generate_partial_diagrams(G)
    directional_diagrams = diagrams.generate_directional_diagrams(G)
    if gen_flux_diagrams:
        flux_diagrams = diagrams.generate_all_flux_diagrams(G)
        flux_diagrams = _flatten_flux_diagrams(flux_diagrams)

    #==========================
    #=== Save the diagrams ====
    #==========================
    # plot and save the input diagram
    plotting.draw_diagrams(G, pos=node_positions, path=save_path, label="input")
    # plot and save the partial diagrams as a panel
    plotting.draw_diagrams(
        partial_diagrams,
        pos=node_positions,
        path=save_path,
        label="partial",
        rows=partial_rows,
    )
    # plot and save the directional diagrams as a panel
    plotting.draw_diagrams(
        directional_diagrams,
        pos=node_positions,
        path=save_path,
        cbt=True,
        label="directional",
        rows=directional_rows,
    )
    if gen_flux_diagrams:
        # plot and save the flux diagrams as a panel
        plotting.draw_diagrams(
            flux_diagrams,
            pos=node_positions,
            path=save_path,
            cbt=True,
            label="flux",
        )
