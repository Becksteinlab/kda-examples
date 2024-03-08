"""
================
Figure 7: Timing
================
"""
import os
from os.path import join, split
import argparse

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker


def fit_powerlaw(diagrams, kda_data):
    time = kda_data[:, 0]

    # make data linear by taking the log10
    log_diagrams = np.log10(diagrams)
    log_time = np.log10(time)

    # get a linear fit of the data (Y = mX + b)
    m, b = np.polyfit(x=log_diagrams, y=log_time, deg=1)
    # for y = ax ** k, log y = k log x + log a
    # translated into linear form, this means
    # m = k, b = log a, X = log x, and Y = log y
    k = m
    a = 10 ** b
    fit_func = a * (diagrams ** k)
    return fit_func, a, k


def get_fit_string(a, k):
    return r"$t_\mathrm{run} = %1.1g D^{%1.3g}$" % (a, k)


def get_dirpar_data(mat_time, kda_time, dirpars):
    unique_dirpars = np.unique(dirpars)
    mat_data = []
    kda_data = []
    for n_dirpars in unique_dirpars:
        mask = dirpars == n_dirpars
        matt = mat_time[mask]
        kdat = kda_time[mask]
        mat_avg = np.mean(matt)
        kda_avg = np.mean(kdat)
        mat_std = np.std(matt)
        kda_std = np.std(kdat)
        mat_data.append([mat_avg, mat_std])
        kda_data.append([kda_avg, kda_std])
    mat_data = np.array(mat_data)
    kda_data = np.array(kda_data)
    return (unique_dirpars, mat_data, kda_data)


def plot_t_over_dirpars(mat_time, kda_time, dirpars, datapath):
    unique_dirpars, mat_data, kda_data = get_dirpar_data(
        mat_time=mat_time,
        kda_time=kda_time,
        dirpars=dirpars,
    )

    # get fit for log-log plot
    fit_func, a, k = fit_powerlaw(diagrams=unique_dirpars, kda_data=kda_data)
    fit_func_str = get_fit_string(a, k)

    fig = plt.figure(figsize=(4, 3), tight_layout=True)
    ax = fig.add_subplot(111)
    ax.set_xscale("log", nonpositive="clip")
    ax.set_yscale("log", nonpositive="clip")
    min_y = 0.5 * np.min(np.abs(mat_data[:, 0] - mat_data[:, 1]))
    max_y = 5 * np.max(np.abs(kda_data[:, 0] + kda_data[:, 1]))
    ax.set_ylim(bottom=min_y, top=max_y)

    # plot fit line
    ax.plot(unique_dirpars, fit_func, color="#7570b3", ls="-", lw=1.5, label=fit_func_str)

    data_arrs = [kda_data, mat_data]
    data_names = ["KDA", "MAT"]
    data_colors = ["#1b9e77", "#d95f02"]
    for data, name, color in zip(data_arrs, data_names, data_colors):
        ax.errorbar(
            unique_dirpars,
            data[:, 0],
            yerr=data[:, 1],
            label=name,
            fmt=".",
            barsabove=False,
            # marker color
            color=color,
            # marker edge color
            mec="black",
            # error bar color
            ecolor="dimgrey",
            # marker edge width
            mew=0.3,
            # error bar line width
            elinewidth=0.5,
            # error bar cap size
            capsize=4,
            # error bar cap thickness
            capthick=0.8,
        )

    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles, labels, loc="upper left", bbox_to_anchor=(-0.02, 1.03))
    ax.set_ylabel("Mean Run Time (s)")
    ax.set_xlabel("Number of Directional Diagrams")

    # set y ticks
    y_ticks = [1e-4, 1e-3, 1e-2, 1e-1, 1e0, 1e1, 1e2]
    ax.set_yticks(y_ticks)
    y_minor = ticker.LogLocator(subs=np.arange(1.0, 10.0) * 0.1, numticks=10)
    ax.yaxis.set_minor_locator(y_minor)
    ax.yaxis.set_minor_formatter(ticker.NullFormatter())

    sns.despine(fig=fig, offset=3)
    fig_savepath = join(datapath, f"time_vs_dirpars.pdf")
    print(f"Saving directional diagrams plot at location: {fig_savepath}")
    fig.savefig(fig_savepath, dpi=500)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("all_data_path", type=str, help="Path to all_data.csv.")
    args = parser.parse_args()
    all_data_path = args.all_data_path

    if os.path.isdir(all_data_path):
        raise Exception(
            "Input path is a directory. Please input a path to all_data.csv."
        )
    print(f"Pulling data from {all_data_path}")

    # get the directory to store generated graphs in
    datapath = split(all_data_path)[0]

    # read in `all_data.csv` and collect required data
    df = pd.read_csv(all_data_path)
    dirpars = df["n_dirpars"]
    mat_time = df["mat time (s)"]
    kda_time = df["kda time (s)"]

    plot_t_over_dirpars(mat_time, kda_time, dirpars, datapath)
