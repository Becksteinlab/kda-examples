import numpy as np
from numpy.testing import assert_allclose
import pandas as pd
from sympy import symbols
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import Patch


def get_var_names_8_state():
    var_names = [
        "k1",
        "k2",
        "k3",
        "k4",
        "k5",
        "k6",
        "k7",
        "k8",
        "k9",
        "k10",
        "k11",
        "k12",
        "k13",
        "k14",
        "k15",
        "k16",
        "k17",
        "k18",
        "k19",
        "k20",
        "k21",
        "k22",
        "k23",
        "k24",
    ]
    return var_names


def data_7b(k_EH_H, k_ED_D, k_EHD_H, k_EHD_D, conf=1):
    """
    Returns the dataset (in tuple form) for Hussey et al. Figure 7B.

    Parameters
    ----------
    k_EH_H : [1, 1e6] --- k13, k75
    k_ED_D : [1e-2, 1e4] --- k43, k65
    k_EHD_H : [1, 1e6] --- k24, k86
    k_EHD_D : [1e-2, 1e4] --- k21, k87

    k_ED_D / k_EHD_D = k_EH_H / k_EHD_H

    Notes
    -----
    Rate mapping (from KDA to Hussey et al):

        k31 <-> k1 <-> H_on
        k13 <-> k2 <-> H_off
        k57 <-> k3 <-> H_on
        k75 <-> k4 <-> H_off
        k42 <-> k5 <-> H_on
        k24 <-> k6 <-> H_off
        k68 <-> k7 <-> H_on
        k86 <-> k8 <-> H_off
        k34 <-> k9 <-> D_on
        k43 <-> k10 <-> D_off
        k56 <-> k11 <-> D_on
        k65 <-> k12 <-> D_off
        k12 <-> k13 <-> D_on
        k21 <-> k14 <-> D_off
        k78 <-> k15 <-> D_on
        k87 <-> k16 <-> D_off
        k71 <-> k17 <-> k_conf_anti
        k17 <-> k18 <-> k_conf_anti
        k53 <-> k19 <-> k_conf_sym
        k35 <-> k20 <-> k_conf_sym
        k64 <-> k21 <-> k_conf_anti
        k46 <-> k22 <-> k_conf_anti
        k82 <-> k23 <-> k_conf_sym
        k28 <-> k24 <-> k_conf_sym
    """
    H_on = 1e10
    D_on = 1e7

    data = (
        H_on,  # 1
        k_EH_H,
        H_on,  # 3
        k_EH_H,
        H_on,  # 5
        k_EHD_H,
        H_on,  # 7
        k_EHD_H,  # in table 3 this is labeled that it has the range of the drug off rates
        D_on,  # 9
        k_ED_D,
        D_on,  # 11
        k_ED_D,
        D_on,  # 13
        k_EHD_D,
        D_on,  # 15
        k_EHD_D,
        conf,  # 17
        conf,
        conf,
        conf,
        conf,
        conf,
        conf,
        conf,  # 24
    )

    return data


def get_7b_dataframe():
    """
    To test how much changing the value of the substrate-off rate
    constants alone could favor symport or antiport phenotypes, we
    varied proton-off rate constants from 1 to 1,000,000 /s and drugoff
    rate constants from 0.01 to 10,000 /s while holding alternating-access
    rate constants uniformly constant at 1, 10, or 100 /s such that RAA = 1

    """
    var_names = get_var_names_8_state()

    # favor antiport initially by using large EHD values first
    EH_H = np.logspace(-2, 8, 201)
    EHD_H = np.logspace(-2, 8, 201)[::-1]
    ED_D = np.logspace(-4, 6, 201)
    EHD_D = np.logspace(-4, 6, 201)[::-1]
    conf_arr = np.array([1, 10, 100, 1000])

    R_off_H = EH_H / EHD_H
    R_off_D = ED_D / EHD_D
    # double check that we are getting identical arrays here
    assert_allclose(R_off_H, R_off_D)

    n_datasets = int(EH_H.size * conf_arr.size)
    data = np.empty((len(var_names), n_datasets))
    i = 0
    for k_conf in conf_arr:
        for (k_EH_H, k_EHD_H, k_ED_D, k_EHD_D) in zip(EH_H, EHD_H, ED_D, EHD_D):
            data[:, i] = data_7b(
                k_EH_H=k_EH_H,
                k_ED_D=k_ED_D,
                k_EHD_H=k_EHD_H,
                k_EHD_D=k_EHD_D,
                conf=k_conf,
            )
            i += 1

    df = pd.DataFrame(columns=var_names, index=np.arange(n_datasets), data=data.T)
    df["figure"] = "7B"
    return df


def sub_dict_7b_values():
    (
        k31,
        k13,
        k57,
        k75,
        k42,
        k24,
        k68,
        k86,
        k34,
        k43,
        k56,
        k65,
        k12,
        k21,
        k78,
        k87,
        k71,
        k17,
        k53,
        k35,
        k64,
        k46,
        k82,
        k28,
    ) = symbols(
        "k31 k13 k57 k75 k42 k24 k68 k86 k34 k43 k56 k65 k12 k21 k78 k87 k71 k17 k53 k35 k64 k46 k82 k28"
    )

    (k_EH_H, k_EHD_H, k_ED_D, k_EHD_D, k_conf) = symbols(
        "k_EH_H k_EHD_H k_ED_D k_EHD_D k_conf"
    )

    H_on = 1e10
    D_on = 1e7
    H_in = 10 ** (-6.5)
    H_out = 10 ** (-7.5)
    c_D = 25e-9

    sub_dict = {
        k31: H_on * H_out,
        k13: k_EH_H,
        k57: H_on * H_in,
        k75: k_EH_H,
        k42: H_on * H_out,
        k24: k_EHD_H,
        k68: H_on * H_in,
        k86: k_EHD_H,
        k34: D_on * c_D,
        k43: k_ED_D,
        k56: D_on * c_D,
        k65: k_ED_D,
        k12: D_on * c_D,
        k21: k_EHD_D,
        k78: D_on * c_D,
        k87: k_EHD_D,
        k71: k_conf,
        k17: k_conf,
        k53: k_conf,
        k35: k_conf,
        k64: k_conf,
        k46: k_conf,
        k82: k_conf,
        k28: k_conf,
    }
    return sub_dict


def rect_label(color):
    rect = Patch(edgecolor="black", facecolor=color, linewidth=0.5, fill=True)
    return rect


def plot_fig_7B(df, colors):

    fig = plt.figure(figsize=(4, 4), tight_layout=True)
    ax = fig.add_subplot(111)  # big subplot
    axH = fig.add_subplot(311)
    axD = fig.add_subplot(312, sharex=axH)

    # Turn off axis lines and ticks of the big subplot
    ax.spines["top"].set_color("none")
    ax.spines["bottom"].set_color("none")
    ax.spines["left"].set_color("none")
    ax.spines["right"].set_color("none")
    ax.tick_params(labelcolor="w", top=False, bottom=False, left=False, right=False)

    column_keys = df.columns
    H_flux_key = column_keys[-2]
    D_flux_key = column_keys[-1]

    R_off_data = df["R_off"].values
    k_AA_data = df["k_AA"].values
    H_flux_data = df[H_flux_key].values
    D_flux_data = df[D_flux_key].values

    unique_k_AA = np.unique(k_AA_data)

    linewidth = 0.6

    axH.axvline(
        x=1,
        ymin=0,
        ymax=1,
        ls="--",
        lw=0.8,
        color="black",
        label=r"$\mathrm{R}_\mathrm{off}$ = 1",
    )
    axD.axvline(
        x=1,
        ymin=0,
        ymax=1,
        ls="--",
        lw=0.8,
        color="black",
        label=r"$\mathrm{R}_\mathrm{off}$ = 1",
    )

    rect_handles = [rect_label(color="black")]
    for i, val in enumerate(unique_k_AA):
        mask = k_AA_data == val
        color = colors[i]
        rect_handles.append(rect_label(color=color))
        label = f"{val:.0e}"
        axH.plot(
            R_off_data[mask],
            H_flux_data[mask],
            ls="-",
            lw=linewidth,
            color=color,
            label=label,
        )
        axD.plot(
            R_off_data[mask],
            D_flux_data[mask],
            ls="-",
            lw=linewidth,
            color=color,
            label=label,
        )

    axH.set_ylabel(r"$\mathrm{H}^{+}$" + " Flux (s$^{-1}$)")
    axD.set_ylabel(r"Drug Flux (s$^{-1}$)")

    _, labels = axH.get_legend_handles_labels()
    ax.legend(
        rect_handles[::-1],
        labels[::-1],
        bbox_to_anchor=(1, 0.5),
        loc="center left",
        title=r"$\mathrm{k}_\mathrm{AA}$",
        handlelength=1,
        handleheight=1,
    )

    max_power = np.max(np.log10(R_off_data))
    xticks = np.logspace(-max_power, max_power, 9)
    logfmt = mpl.ticker.LogFormatterSciNotation(base=10.0, labelOnlyBase=True)
    for ax in [axH, axD]:
        ax.grid()
        ax.set_xscale("log")
        ax.set_xticks([])
        ax.set_xticklabels([])
        ax.set_xticks(xticks)
        ax.set_xticklabels(xticks)
        ax.xaxis.set_major_formatter(logfmt)

    axH.tick_params(axis="both", which="both", labelsize=7, labelbottom=False)
    axD.tick_params(axis="both", which="both", labelsize=7)

    return fig, [axH, axD]
