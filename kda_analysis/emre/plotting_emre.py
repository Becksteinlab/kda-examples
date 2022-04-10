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


def is_odd(x):
    return x % 2 != 0


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


def plot_cycle_flux(RAA, cycles, fluxes, flux_tot, species, raw=False):
    colors = get_colors()
    fig = plt.figure(figsize=(13, 7), tight_layout=True)
    ax = fig.add_subplot(111)

    if species == "H":
        label1 = "Proton"
        label2 = "J_{H}"
        colour = "#A02020"
    elif species == "D":
        label1 = "Drug"
        label2 = "J_{D}"
        colour = "#6BE35D"
    if raw == False:
        i = 1
        for cycle, cf, col in zip(cycles, fluxes, colors):
            new_cycle = list(np.array(cycle) + 1)
            if is_odd(i) == True:
                style = "-"
            else:
                style = "--"
            ax.semilogx(
                RAA,
                100 * np.abs(cf) / np.abs(flux_tot),
                ls=style,
                lw=2,
                color=col,
                label="{}".format(new_cycle),
            )
            i += 1
        ax.set_ylim(0, 100)
        ax.set_title("EmrE {} Flux Contributions by Cycle".format(label1))
        ax.set_ylabel(r"Contribution (%)")

    elif raw == True:
        ax.semilogx(
            RAA, flux_tot, "-", lw=2, color=colour, label=r"${}$".format(label2)
        )
        i = 1
        for cycle, cf, col in zip(cycles, fluxes, colors):
            new_cycle = list(np.array(cycle) + 1)
            if is_odd(i) == True:
                style = "-"
            else:
                style = "--"
            ax.semilogx(
                RAA, cf, ls=style, lw=2, color=col, label="{}".format(new_cycle)
            )
            i += 1
        ax.set_title("EmrE {} Flux by Cycle".format(label1))
        ax.set_ylabel(r"Flux ($s^{-1}$)")

    ax.axhline(y=0, xmin=0, xmax=1, ls="-", color="black")
    ax.axvline(x=1e0, ymin=0, ymax=1, ls="-", color="black")
    ax.set_xlabel(r"$R_{AA}$")
    plt.legend(bbox_to_anchor=(1, 1), loc="upper left", ncol=1)
    ax.grid(True)


def plot_transition_flux(RAA, edges, fluxes, species):
    colors = get_colors()
    fig = plt.figure(figsize=(12, 7), tight_layout=True)
    ax = fig.add_subplot(111)

    if species == "H":
        label1 = "Proton"
        label2 = "J_{H}"
        colour = "#A02020"
    elif species == "D":
        label1 = "Drug"
        label2 = "J_{D}"
        colour = "#6BE35D"

    i = 1
    for edge, tf, col in zip(edges, fluxes, colors):
        new_edge = [edge[0], edge[1]]
        edge_tuple = tuple(np.array(new_edge) + 1)
        if is_odd(i) == True:
            style = "-"
        else:
            style = "--"
        ax.semilogx(RAA, tf, ls=style, lw=2, color=col, label="J_{}".format(edge_tuple))
        i += 1
    ax.set_title("EmrE {} Flux by Transition".format(label1))
    ax.set_ylabel(r"Flux ($s^{-1}$)")

    ax.axhline(y=0, xmin=0, xmax=1, ls="-", color="black")
    ax.axvline(x=1e0, ymin=0, ymax=1, ls="-", color="black")
    ax.set_xlabel(r"$R_{AA}$")
    plt.legend(bbox_to_anchor=(1, 1), loc="upper left", ncol=1)
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
