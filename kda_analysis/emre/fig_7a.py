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


def data_7a(k_AA_anti=1, k_AA_sym=1):
    """
    Returns the dataset (in tuple form) for Hussey et al. Figure 7A. Since
    they used a range of conformational change rates (uniformly), this takes
    in a parameter that allows that value to be changed.

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
    H_off = 1e3
    D_on = 1e7
    D_off = 1
    data = (
        H_on,
        H_off,
        H_on,
        H_off,
        H_on,
        H_off,
        H_on,
        H_off,
        D_on,
        D_off,
        D_on,
        D_off,
        D_on,
        D_off,
        D_on,
        D_off,
        k_AA_anti,
        k_AA_anti,
        k_AA_sym,
        k_AA_sym,
        k_AA_anti,
        k_AA_anti,
        k_AA_sym,
        k_AA_sym,
    )

    return data


def get_7a_dataframe():
    var_names = get_var_names_8_state()

    # according to https://github.com/granthussey/EmrE-Modeling/blob/master/figure7a.m
    # the rates were generated using a logarithmic scale:
    # round(logspace(0,2,100),1);  % Logspace vector from 1 to 100 rounded to the nearest 10ths place.

    # don't round values since some of the array elements end up
    # being the same value, and generate 101 elements so we hit the
    # integers 1, 10, and 100 precisely
    n_datasets = 201
    # favor antiport initially
    k_AA_anti = np.logspace(-3, 5, n_datasets)[::-1]
    k_AA_sym = np.logspace(-3, 5, n_datasets)
    # generate an empty array to store the data set
    data = np.empty((len(var_names), n_datasets))
    for i, (k_anti, k_sym) in enumerate(zip(k_AA_anti, k_AA_sym)):
        data[:, i] = data_7a(k_AA_anti=k_anti, k_AA_sym=k_sym)
    df = pd.DataFrame(columns=var_names, index=np.arange(n_datasets), data=data.T)
    df["figure"] = "7A"
    return df


def sub_dict_7a_symbols():

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

    (H_on, H_off, D_on, D_off, H_in, H_out, D_in, D_out, k_AA_anti, k_AA_sym) = symbols(
        "H_on H_off D_on D_off H_in H_out D_in D_out k_AA_anti k_AA_sym"
    )

    sub_dict = {
        k31: H_on * H_out,
        k13: H_off,
        k57: H_on * H_in,
        k75: H_off,
        k42: H_on * H_out,
        k24: H_off,
        k68: H_on * H_in,
        k86: D_off,
        k34: D_on * D_out,
        k43: D_off,
        k56: D_on * D_in,
        k65: D_off,
        k12: D_on * D_out,
        k21: D_off,
        k78: D_on * D_in,
        k87: D_off,
        k71: k_AA_anti,
        k17: k_AA_anti,
        k53: k_AA_sym,
        k35: k_AA_sym,
        k64: k_AA_anti,
        k46: k_AA_anti,
        k82: k_AA_sym,
        k28: k_AA_sym,
    }
    return sub_dict


def sub_dict_7a_values():
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

    (k_AA_anti, k_AA_sym) = symbols("k_AA_anti k_AA_sym")

    H_on = 1e10
    H_off = 1e3
    D_on = 1e7
    D_off = 1
    H_in = 10 ** (-6.5)
    H_out = 10 ** (-7.5)
    c_D = 25e-9

    sub_dict = {
        k31: H_on * H_out,
        k13: H_off,
        k57: H_on * H_in,
        k75: H_off,
        k42: H_on * H_out,
        k24: H_off,
        k68: H_on * H_in,
        k86: D_off,
        k34: D_on * c_D,
        k43: D_off,
        k56: D_on * c_D,
        k65: D_off,
        k12: D_on * c_D,
        k21: D_off,
        k78: D_on * c_D,
        k87: D_off,
        k71: k_AA_anti,
        k17: k_AA_anti,
        k53: k_AA_sym,
        k35: k_AA_sym,
        k64: k_AA_anti,
        k46: k_AA_anti,
        k82: k_AA_sym,
        k28: k_AA_sym,
    }
    return sub_dict


def plot_fig_7A(df):
    fig = plt.figure(figsize=(4, 3), tight_layout=True)
    ax = fig.add_subplot(111)

    column_keys = df.columns
    H_flux_key = column_keys[-2]
    D_flux_key = column_keys[-1]

    H_colour = "#A02020"
    D_colour = "#6BE35D"
    ax.semilogx(
        df["R_AA"].values,
        df[H_flux_key].values,
        ls="-",
        color=H_colour,
        label=r"$\mathrm{J}_\mathrm{H}$",
    )
    ax.semilogx(
        df["R_AA"].values,
        df[D_flux_key].values,
        ls="-",
        color=D_colour,
        label=r"$\mathrm{J}_\mathrm{D}$",
    )

    ax.axvline(
        x=1,
        ymin=0,
        ymax=1,
        ls="--",
        color="black",
        label=r"$\mathrm{R}_\mathrm{AA}$ = 1",
    )
    plt.legend(bbox_to_anchor=(1, 1), loc="upper left", ncol=1)
    return fig, [ax]
