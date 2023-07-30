"""
Module for building figures from https://doi.org/10.1085/jgp.201912437
"""

import sys
import argparse

import pandas as pd
from sympy import lambdify
import matplotlib.pyplot as plt

from thermo_consistency import run_thermo_consistency_checks
import fig_utils


def main(fig_list, check_thermo_con):
    K, G = fig_utils.get_K_and_G()
    H_cos_dict, D_cos_dict = fig_utils.get_H_D_cos_dicts(G=G, K=K)

    data_path = "data/all_figure_data.csv"
    df_all = pd.read_csv(data_path)

    for fig_key in fig_list:

        print("=" * 45)
        print(f"Generating Fig. {fig_key}")
        print("=" * 45)
        sys.stdout.flush()

        mask = df_all["figure"] == fig_key
        df_fig = df_all[mask]

        if check_thermo_con:
            # NOTE: data_path set to None means it will not update
            # figure_data.csv with the thermo consistency data.
            # To do that, run:
            #    `python check_thermo_consistency.py`
            # on the command line
            run_thermo_consistency_checks(df=df_fig, data_path=None)

        sub_dict = fig_utils.get_sub_dict(fig_key=fig_key)
        rate_names = fig_utils.get_rate_names(fig_key=fig_key)

        H_sympy_func_path = f"fig_functions/fig_{fig_key}_H_sympy.txt"
        D_sympy_func_path = f"fig_functions/fig_{fig_key}_D_sympy.txt"

        print("--> Generating sympy function for H")
        sys.stdout.flush()
        H_op_flux_f = fig_utils.get_sympy_func(
            G, H_cos_dict, sub_dict, filepath=H_sympy_func_path
        )
        print("--> Generating sympy function for D")
        sys.stdout.flush()
        D_op_flux_f = fig_utils.get_sympy_func(
            G, D_cos_dict, sub_dict, filepath=D_sympy_func_path
        )

        # now that we have the simplified SymPy expressions, we can
        # convert them into Python lambda functions
        H_op_flux_f = lambdify(rate_names, H_op_flux_f, "numpy")
        D_op_flux_f = lambdify(rate_names, D_op_flux_f, "numpy")

        plot_df = fig_utils.get_plot_data(
            df=df_fig, fH=H_op_flux_f, fD=D_op_flux_f, fig_key=fig_key
        )
        plot_data_path = f"./data/fig_{fig_key}_data.csv"
        plot_df.to_csv(plot_data_path, index=False)


        fig = fig_utils.plot_data(df=plot_df, fig_key=fig_key)
        save_path = f"plots/figures/fig_{fig_key}"
        print(f"--> Saving Fig. {fig_key} at location: {save_path}...")
        for _ext in ("png", "pdf", "svg"):
            fig.savefig(f"{save_path}.{_ext}", dpi=300)
        plt.close()

        supported_figs = ["7A", "7B"]
        if fig_key in supported_figs:
            data_dict = fig_utils.plot_flux_graphs(df=df_fig, fig_key=fig_key)
            for i, (flux_graph, filename) in enumerate(
                zip(data_dict["figures"], data_dict["filenames"])
            ):
                save_path = f"plots/flux_graphs/{filename}"
                print(f"--> Saving Fig. {fig_key} at location: {save_path}...")
                for _ext in ("png", "pdf", "svg"):
                    flux_graph.savefig(f"{save_path}.{_ext}", dpi=300)
                plt.close()


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-f",
        "--figures",
        nargs="*",
        help=f"Specify figures to plot ('7A', '7B', '7D' or '9'). \n "
        f"If argument not used, will generate all.",
    )

    parser.add_argument(
        "-tc",
        "--check_thermo_con",
        help="turn thermodynamic consistency checking on",
        action="store_true",
    )
    args = parser.parse_args()
    fig_list = args.figures
    check_thermo_con = args.check_thermo_con
    accepted_keys = ["7A", "7B", "7D", "9"]

    if fig_list is None:
        print("No figures specified, running all.")
        fig_list = accepted_keys
    elif len(fig_list) == 0:
        raise ValueError(
            "When using '-f', please specify a figure (i.e. '7A', '7B', etc.)"
        )
    else:
        for fig_key in fig_list:
            if not fig_key in accepted_keys:
                raise NotImplementedError(
                    f"Input figure {fig_key} is not supported. \n"
                    f"Please enter from the following: {accepted_keys}"
                )

    sys.stdout.flush()
    main(fig_list=fig_list, check_thermo_con=check_thermo_con)
