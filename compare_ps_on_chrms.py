import numpy as np
import matplotlib.pyplot as plt
import logging

logging.basicConfig(level=logging.INFO)
import pandas as pd
from fit_functions import plot_ps, fit_linear_regression, fit_ps_log_bins
from plot_functions import multiplots, multiplot_with_subplots
from Expected_calculator import dump
import datetime


def process(data):
    # normalize probabilities to 1
    for chr in data:
        # last expected values on chrms are similar
        # let's clip expected values at this point
        t = np.subtract(data[chr][1:], data[chr][:-1])
        t = np.where(t == 0)[0]
        if len(t) != 0:
            t = max(t)
            assert t > 3 * len(data) / 4
            data[chr] = data[chr][:t]

        # now normalize the sum to 1
        data[chr] = data[chr] / np.sum(data[chr])
    return data


def row2color(row):
    colors = {"Anopheles": "springgreen",
              "Drosophila": "springgreen",
              #            "culex": "limegreen",
              "culex": "springgreen",
              "Aedes": "springgreen",
              "Polypedium": "springgreen",
              "chick": "springgreen",
              "mammals": "springgreen",
              "chick_Dekker": "springgreen"
              }

    special_styles = {
        "Mature erythrocytes": {"linestyle": "-", "color": "blue", "linewidth": 2},
        "CIE": {"marker": "*", "linestyle": ":"},

        "Aedes": {"color": "orange", "linewidth": 1},
        "Culex": {"color": "lightseagreen", "linewidth": 1},
        "HCT116 RAD21-": {"color": "blue", "linewidth": 2},
        "HAP1 WAPL1-": {"color": "red", "linewidth": 2},
        #            "HAP1": {"color":  "red", "linestyle":"--", "linewidth":0.5},
        #            "LiverTAM": {"color": "yellow", "linestyle": "--", "linewidth": 1},
        "Liver Nipbl-": {"color": "red", "linestyle": "--", "linewidth": 1},
        #            "BonevCN":{"color":"purple"},
        #            "DekkerCapH-": {"color":"red","marker" : "*", "linestyle":"--"},
        #            "DekkerCAPHControl": {"color": "red", "marker" : "*"},
        #            "DekkerCapH2-": {"color":"salmon", "marker":"^", "linestyle" : "--"},
        #            "DekkerCAPH2Control": {"color": "salmon", "marker":"^"},
        #            "DekkerSMC2-":{"color": "red"},
        #            "DekkerSMC2Control": {"color": "red", "linestyle":"--"},
        "Chick Prometaphase": {"color": "yellow", "linewidth": 2},
        "Chick SMC2-": {"linestyle": "-", "color": "red", "linewidth": 1},
        "Chick CapH-": {"linestyle": "-", "color": "red", "linewidth": 1},
        "Chick CapH2-": {"linestyle": "-", "color": "red", "linewidth": 1},
        #            "DekkerCapH-":{"linestyle":"--", "color":"red", "linewidth":1},
        #            "DekkerCapH2-":{"linestyle":":", "color":"red", "linewidth":1},
        "Dmel CAPH2-": {"linestyle": "-", "color": "red", "linewidth": 2},
        "Dmel RAD21-": {"linestyle": "-", "color": "blue", "linewidth": 2},
        #            "SextonDrosophila":{"linestyle": "-", "color":"black"},
        #            "Kc167rowley": {"linestyle": "-", "color": "blue"},
        #            "S2": {"linestyle": "-", "color": "red"},
        #            "S2HeatShock": {"linestyle": "--", "color": "red"}
        #        "Dvir": {"color": "white"},
        #        "Dmel": {"color": "blue"}
        #        "Dbus": {"color": "white"},
    }

    result = {"linewidth": 0.5}
    result["color"] = colors[row.subtaxon]
    if row["name"] in special_styles:
        for k, v in special_styles[row["name"]].items():
            result[k] = v

    return result


# this file contains list of data and defines path to .hic files
dataset = "datasets.csv"

# data is stored in .hic format, we need juicer tools to work with it
# jucer tools could be downloaded from https://github.com/aidenlab/juicer/wiki/Download
# after download, specify path to juicer_tools jar file
juicer_tools = "juicer/juicer_tools_1.19.02.jar"

# which resolution to use?
resolution = 10000

# stop analysis after this number of datasets. Used for debug only
report = np.inf

datasets = pd.read_csv(dataset, sep="\t",
                       comment="#")

# each item of this dict describes subset of datasets for analysis,
# e.g. one could subset datasets by taxon
analysis = {
    "Anopheles": datasets.query(
        "(name in ['An. col','An. mer (adult)','An. mer (embryo)','An. ste','An. alb','An. atr'])"),
    #    "test": datasets.query("(name in ['Acol','Amer'])")
    #    "Other_insects": datasets.query("(subtaxon=='Drosophila' or subtaxon=='culex' or name=='Aedes')")
    #    "Other_insects_for_multiplot": datasets.query("(subtaxon=='Drosophila' or subtaxon=='culex' or name=='Aedes')" + \
    #            "and not (name in ['Dmel embryo 3-4h','Dmel embryo nc12','Dmel embryo nc13','Dmel embryo nc1-4','Dmel embryo nc1-4 (mitotic)','Salivary glands(Polytene)'])")
    #    "mammals": datasets.query("(subtaxon=='mammals')")
    #    "chicken": datasets.query("(subtaxon=='chick')")
    #    "Nipbl": datasets.query("(name in ['LiverWT','LiverTAM','LiverNipbl'])")
    #    "Aedes":  datasets.query("name=='Aedes'")
    #    "mammals_test": datasets.query("(name=='BonevNPC')")
    #    "all_maps_from_Gibcus_et_al": datasets.query("subtaxon=='chick_Dekker'")
}

# draw all graphs on one plot (=Trus)
# or draw multiple subplots, one subplot for each dataset (=False)
multiplot = True  # draw all graphs on one plot or draw multiple subplots

# functions are defined in fit_functions.pu
# I mostly use Slope, which computes slopoe, and Ps, which draws P(s)

def chr2color(chr):
    if chr == "X":
        return {"color":"red"}
    return {"color":"black"}

norm="None" # choose KR or None

# for func in ["Ps_log"]:
for func in ["Ps"]:
    # for func in ["Slope"]:
    # for func in funcs:
    for suffix, species in analysis.items():
        # example: suffix = Anopheles, items is pd table with anopheles datasets
        plots = {}  # here we will store computed data for plot
        for ind in range(len(species)):
            row = species.iloc[ind]
            logging.info(row["name"])
            hic = row.link
            logging.info("Starting dump")

            # dump expected vector from juicebox using juicer_tools soft
            # see Expected_calculator.py for details
            data = dump(hic, juicer_tools, resolution, excludechrms=[], norm=norm)
            for chr in data:
                logging.info("Fitting...")

                # fit is the core of whole module. It fits linear regression in log/log coordinates,
                # or returns P(s) values
                # see fit_functions code
                if func == "Ps":
                    X, Y = plot_ps({chr: data[chr]}, resolution=resolution)
                else:
                    raise
                logging.info("Plotting...")
                # plots is a dictionary, each item (type=list) describes one plot
                # description should contain:
                # 1.) pd.Dataframe with X and Y values - this is main data, X and Y points
                # 2.) optional arguments, such as line color, line width and etc.,
                # 3.) legend (=plot name)
                # 4.) sample genotype (currently not used)
                plots[chr] = [pd.DataFrame({"X": X, "Y": Y}), chr2color(chr),
                              row["name"]+chr, row["Genotype"]]
                logging.info("Done!")
                if ind >= report and ind % report == 0:
                    # plt.show()
                    break
            # now when all chrms processed, draw a plot
            # see multiplots and multiplot_with_subplots code in plot_functions.py
            if multiplot:
                multiplots(plots, shadow=(func == "Slope"), average=(func == "Slope"))
                # multiplots(plots, shadow=False, average=False)
                if func == "Slope":
                    plt.gca().set_ylabel("Slope")
                elif func == "Ps":
                    plt.gca().set_ylabel("Log(Normed Contact probability)")
                plt.gca().set_xlabel("Genomic distance")
            else:
                multiplot_with_subplots(plots, xlabel="Genomic distance", y_label="Slope")
            # plt.tight_layout()
            fig_path = "results/" + str(datetime.datetime.today().date()) + \
                       "result_" + suffix + "_" + row["name"] + "_norm_" + norm + "_" + func + "_" + str(multiplot) + ".png"
            print("Saving figure " + fig_path)
            plt.savefig(fig_path, dpi=500)
            plt.clf()