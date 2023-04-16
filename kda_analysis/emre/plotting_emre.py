import numpy as np
import matplotlib
import matplotlib.pyplot as plt

def get_colors():
    color_list = [
        "#e6194B",
        "#3cb44b",
        "#ffe119",
        "#4363d8",
        "#f58231",
        "#911eb4",
        "#42d4f4",
        "#f032e6",
        "#bfef45",
        "#fabed4",
        "#469990",
        "#9A6324",
        "#800000",
        "#808000",
        "#000075",
        "#a9a9a9",
        "#000000",
    ]
    return color_list


def get_node_positions():
    node_pos = {
        0: [2, 2],
        1: [1, 1],
        2: [-2, 2],
        3: [-1, 1],
        4: [-2, -2],
        5: [-1, -1],
        6: [2, -2],
        7: [1, -1],
    }
    return node_pos


def get_linestyle(x):
    if x % 2 != 0:
        return "-"
    else:
        return "--"


def get_linewidth(x):
    if x % 2 != 0:
        return 1.25
    else:
        return 2.5


def plot_net_cycle_fluxes(R, fluxes, flux_tot, cycles, norm=True):

    fig = plt.figure(figsize=(5, 4), tight_layout=True)
    ax = fig.add_subplot(111)

    # if any of the net cycle flux values across the range 
    # are above this threshold they will be plotted
    threshold = 1e-3
    color_list = get_colors()
    color_idx = 0

    if norm:
        for i, (cycle, flux) in enumerate(zip(cycles, fluxes)):
            if np.any(flux > threshold):
                normalized_flux = 100 * np.abs(flux) / np.abs(flux_tot)
                ax.semilogx(
                    R,
                    normalized_flux,
                    color=color_list[color_idx],
                    ls=get_linestyle(color_idx),
                    lw=get_linewidth(color_idx),
                    label=f"Cycle {i+1}",
                )
                color_idx +=1
            ax.set_ylim(0, 100)
            ax.set_title("EmrE Net Cycle Flux % Contribution")
            ax.set_ylabel("Contribution (%)")

    else:
        ax.semilogx(x, flux_tot, "-", lw=2, color="black", label="Total")
        
        for i, (cycle, flux) in enumerate(zip(cycles, fluxes)):
            if np.any(flux > threshold):
                ax.semilogx(
                    R, 
                    flux,
                    color=color_list[color_idx], 
                    ls=get_linestyle(color_idx), 
                    lw=get_linewidth(color_idx), 
                    label=f"Cycle {i+1}",
                )
                color_idx +=1
        ax.set_title("EmrE Net Cycle Fluxes")
        ax.set_ylabel(r"Flux ($s^{-1}$)")

    ax.axvline(x=1e0, ymin=0, ymax=1, ls="--", label=r"$R_{\mathrm{AA}}=1$", color="black")
    ax.set_xlabel(r"$R_{\mathrm{AA}}$")
    plt.legend(bbox_to_anchor=(1, 1), loc="upper left", ncol=1)
    ax.grid(True)
    plt.close()
    return fig


def plot_transition_flux(RAA, edges, fluxes):
    colors = get_colors()
    fig = plt.figure(figsize=(8, 5), tight_layout=True)
    ax = fig.add_subplot(111)

    threshold = 1e-3

    color_idx = 0
    for edge, flux, col in zip(edges, fluxes, colors):
        if np.any(flux > threshold):
            new_edge = [edge[0], edge[1]]
            edge_tuple = tuple(np.array(new_edge) + 1)
            ax.semilogx(
                RAA, 
                flux, 
                ls=style, 
                lw=2, 
                color=col, 
                label=f"J_{edge_tuple}",
            )
            color_idx += 1
    ax.set_title("EmrE Net Transition Fluxes")
    ax.set_ylabel(r"Flux ($s^{-1}$)")

    ax.axhline(y=0, xmin=0, xmax=1, ls="-", color="black")
    ax.axvline(x=1e0, ymin=0, ymax=1, ls="-", color="black")
    ax.set_xlabel(r"$R_{\mathrm{AA}}$")
    plt.legend(bbox_to_anchor=(1, 1), loc="upper left", ncol=1)
    ax.grid(True)



def plot_gradients(RAA, Tr, Hr):
    fig = plt.figure(figsize=(8, 6), tight_layout=True)
    ax = fig.add_subplot(111)
    ax.loglog(
        RAA,
        Tr,
        ls="-",
        lw=2,
        color="#6BE35D",
        label=r"$T_r = \frac{[D]_{in}}{[D]_{out}}$",
    )
    ax.loglog(
        RAA,
        Hr,
        ls="-",
        lw=2,
        color="#A02020",
        label=r"$H_r = \frac{[H]_{in}}{[H]_{out}}$",
    )
    ax.set_title("EmrE Concentration Gradients")
    ax.set_ylabel(r"Gradient")
    ax.set_xlabel(r"$R_{AA}$")
    # ax.axhline(y=1e0, xmin=0, xmax=1, ls='-', color='black')
    # ax.axvline(x=1e0, ymin=0, ymax=1, ls='-', color='black')
    plt.legend(loc="best")
    ax.grid(True)


def plot_3D_flux(X, Y, Z, vmin, vmax, azimuth=-50, elev=40):
    fig = plt.figure(figsize=(9, 6), tight_layout=True)
    ax = fig.add_subplot(111, projection="3d")
    cmap = plt.cm.bwr
    norm = matplotlib.colors.Normalize(vmin, vmax)
    offset = 0
    zlabel = r"Flux ($s^{-1}$)"
    ax.set_xlabel(r"log($R_{AA}$)")
    ax.set_ylabel(r"log($T_{r}$)")
    ax.set_zlim(Z.min(), Z.max())
    ax.view_init(elev=elev, azim=azimuth)
    ax.plot_wireframe(
        np.log10(X),
        np.log10(Y),
        Z,
        rstride=10,
        cstride=10,
        linewidth=0.8,
        color="black",
    )
    surf = ax.plot_surface(
        np.log10(X), np.log10(Y), Z, cmap=cmap, norm=norm, alpha=0.65
    )
    cset = ax.contourf(
        np.log10(X), np.log10(Y), Z, 10, zdir="z", offset=Z.min(), cmap=cmap, norm=norm
    )
    cb = fig.colorbar(surf, cmap=cmap, norm=norm, shrink=0.5, aspect=7)
    cb.set_label(zlabel)
