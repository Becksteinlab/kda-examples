import numpy as np
import pandas as pd
from sympy import symbols
import matplotlib.pyplot as plt


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


def data_9():
    """
    Returns the dataset (in tuple form) for Hussey et al. Figure 9.

    Parameters
    ----------
    k_AA_anti : int
        The conformational change rate used for (symport) rates
        k19, k20, k23, k24. Should be in the range [1, 100].
    k_AA_sym : int
        The conformational change rate used for (symport) rates
        k17, k18, k21, k22. Should be in the range [1, 100].
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
    H_off_anti = 63.1
    H_off_sym = 1584

    D_on = 1e7
    D_off_anti = 0.631
    D_off_sym = 10

    k_AA_anti_1 = 100
    k_AA_anti_2 = 7.3

    k_AA_sym_1 = 40
    k_AA_sym_2 = 8.9

    data = np.array(
        [
            H_on,  # k31 <-> k1 <-> H_on
            H_off_anti,  # k13 <-> k2 <-> H_off
            H_on,  # k57 <-> k3 <-> H_on
            H_off_anti,  # k75 <-> k4 <-> H_off
            H_on,  # k42 <-> k5 <-> H_on
            H_off_sym,  # k24 <-> k6 <-> H_off
            H_on,  # k68 <-> k7 <-> H_on
            H_off_sym,  # k86 <-> k8 <-> H_off
            D_on,  # k34 <-> k9 <-> D_on
            D_off_anti,  # k43 <-> k10 <-> D_off
            D_on,  # k56 <-> k11 <-> D_on
            D_off_anti,  # k65 <-> k12 <-> D_off
            D_on,  # k12 <-> k13 <-> D_on
            D_off_sym,  # k21 <-> k14 <-> D_off
            D_on,  # k78 <-> k15 <-> D_on
            D_off_sym,  # k87 <-> k16 <-> D_off
            k_AA_anti_1,  # k71 <-> k17 <-> k_conf_anti
            k_AA_anti_1,  # k17 <-> k18 <-> k_conf_anti
            k_AA_sym_1,  # k53 <-> k19 <-> k_conf_sym
            k_AA_sym_1,  # k35 <-> k20 <-> k_conf_sym
            k_AA_anti_2,  # k64 <-> k21 <-> k_conf_anti
            k_AA_anti_2,  # k46 <-> k22 <-> k_conf_anti
            k_AA_sym_2,  # k82 <-> k23 <-> k_conf_sym
            k_AA_sym_2,  # k28 <-> k24 <-> k_conf_sym
        ]
    )

    return data


def get_9_dataframe():
    # the rate constants don't change here, so we just need to call the
    # `data_9` function and convert it to a pandas dataframe
    var_names = get_var_names_8_state()

    data = data_9()
    df = pd.DataFrame(columns=[0], index=var_names, data=data)
    df = df.T
    df["figure"] = "9"
    return df


def sub_dict_9_symbols():

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

    (
        H_int,
        H_ext,
        D_int,
        D_ext,
        H_on,
        H_off_anti,
        H_off_sym,
        D_on,
        D_off_anti,
        D_off_sym,
        k_AA_anti_1,
        k_AA_anti_2,
        k_AA_sym_1,
        k_AA_sym_2,
    ) = symbols(
        "H_int H_ext D_int D_ext H_on H_off_anti H_off_sym D_on D_off_anti D_off_sym k_AA_anti_1 k_AA_anti_2 k_AA_sym_1 k_AA_sym_2"
    )

    sub_dict = {
        k31: H_on * H_ext,
        k13: H_off_anti,
        k57: H_on * H_int,
        k75: H_off_anti,
        k42: H_on * H_ext,
        k24: H_off_sym,
        k68: H_on * H_int,
        k86: H_off_sym,
        k34: D_on * D_ext,
        k43: D_off_anti,
        k56: D_on * D_int,
        k65: D_off_anti,
        k12: D_on * D_ext,
        k21: D_off_sym,
        k78: D_on * D_int,
        k87: D_off_sym,
        k71: k_AA_anti_1,
        k17: k_AA_anti_1,
        k53: k_AA_sym_1,
        k35: k_AA_sym_1,
        k64: k_AA_anti_2,
        k46: k_AA_anti_2,
        k82: k_AA_sym_2,
        k28: k_AA_sym_2,
    }
    return sub_dict


def sub_dict_9_values():

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

    # the only value they vary here is the external pH value (-log10([H]/[H0]))
    H_ext = symbols("H_ext")

    H_int = 10 ** (-7.4)
    D_int = 25e-9
    D_ext = 25e-9
    H_on = 1e10
    H_off_anti = 63.1
    H_off_sym = 1584
    D_on = 1e7
    D_off_anti = 0.631
    D_off_sym = 10
    k_AA_anti_1 = 100
    k_AA_anti_2 = 7.3
    k_AA_sym_1 = 40
    k_AA_sym_2 = 8.9

    sub_dict = {
        k31: H_on * H_ext,
        k13: H_off_anti,
        k57: H_on * H_int,
        k75: H_off_anti,
        k42: H_on * H_ext,
        k24: H_off_sym,
        k68: H_on * H_int,
        k86: H_off_sym,
        k34: D_on * D_ext,
        k43: D_off_anti,
        k56: D_on * D_int,
        k65: D_off_anti,
        k12: D_on * D_ext,
        k21: D_off_sym,
        k78: D_on * D_int,
        k87: D_off_sym,
        k71: k_AA_anti_1,
        k17: k_AA_anti_1,
        k53: k_AA_sym_1,
        k35: k_AA_sym_1,
        k64: k_AA_anti_2,
        k46: k_AA_anti_2,
        k82: k_AA_sym_2,
        k28: k_AA_sym_2,
    }
    return sub_dict


def plot_fig_9(df):
    fig = plt.figure(figsize=(4, 4), tight_layout=True)
    ax = fig.add_subplot(211)

    column_keys = df.columns
    pH_key = column_keys[0]
    H_flux_key = column_keys[-2]
    D_flux_key = column_keys[-1]

    H_colour = "#A02020"
    D_colour = "#6BE35D"
    ax.semilogx(
        df[pH_key].values,
        df[H_flux_key].values,
        ls="-",
        color=H_colour,
        label=r"$J_\mathrm{H}$",
    )
    ax.semilogx(
        df[pH_key].values,
        df[D_flux_key].values,
        ls="-",
        color=D_colour,
        label=r"$J_\mathrm{D}$",
    )

    x_ticks = np.linspace(5.4, 9.4, 5)
    ax.set_xticks(x_ticks, minor=False)
    ax.set_xticks([], minor=True)
    ax.set_xticklabels(x_ticks)

    ax.axvline(x=7.4, ymin=0, ymax=1, ls="--", color="black")
    return fig, [ax]
