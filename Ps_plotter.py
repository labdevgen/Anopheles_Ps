import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import subprocess
import straw
import logging
logging.basicConfig(level = logging.INFO)
from hashlib import md5
import os
import pandas as pd
from sklearn.linear_model import LinearRegression
import pickle

# dump data to .expected file or reload from cash
def dump(file, juicerpath, resolution, minchrsize = 20000000, excludechrms = ("X","chrX",
                                                                              "Y","chrY",
                                                                              "Z","chrZ",
                                                                              "W","chrW")):
    description = (file+str(resolution)+str(minchrsize)+str(excludechrms)).encode("utf-8")

    hashfile = os.path.join("data/all_chr_dumps",
                            md5(description).hexdigest())
    if os.path.isfile(hashfile):
        return pickle.load(open(hashfile,"rb"))

    # get chrm names first
    strawObj = straw.straw(file)
    chrsizes = strawObj.chromDotSizes.data
    # example:
    # {'ALL': (0, 231617), 'X': (1, 26840812), '2R': (2, 61076761), '2L': (3, 48179558), '3R': (4, 53444961),
    # '3L': (5, 42074947)}

    if "ALL" in chrsizes:
        del chrsizes["ALL"]

    # get valid chrms
    valid_chrms = []

    for chr in chrsizes:
        if chrsizes[chr][1] > minchrsize and \
                not chr in excludechrms and \
                not chr.lower() in excludechrms and \
                not chr.upper() in excludechrms:
            valid_chrms.append(chr)

    assert len(valid_chrms) > 0

    # maxlenid = sorted(chrsizes.keys(),
    #                   key=lambda x: chrsizes[x][1],
    #                   reverse=True)[0]
    # valid_chrms = [maxlenid]
    # logging.info("Using chr "+str(maxlenid)+" for file "+file)

    # dump data for all valid chrms
    data = {}
    for chr in valid_chrms:
        description_chr = (file+str(resolution)+str(minchrsize)+str(excludechrms)+chr).encode("utf-8")
        chromosome_hash_file = os.path.join("data/all_chr_dumps/",
                            md5(description_chr).hexdigest())
        if not os.path.isfile(chromosome_hash_file):
            command = ["java","-jar",juicerpath,"dump","expected",
                        "KR",file,chr,"BP",str(resolution),chromosome_hash_file]
            logging.info(" ".join(command))
            p = subprocess.run(" ".join(command),shell=True,
                                    stdout = subprocess.PIPE,
                                    stderr=subprocess.PIPE,
                                    check=True)
        chr_data = np.loadtxt(chromosome_hash_file)
        data[chr] = chr_data

    pickle.dump(data, open(hashfile,"wb"))
    return data

def process(data):
    # normalize probabilities to 1
    for chr in data:
        t = np.subtract(data[chr][1:], data[chr][:-1])
        t = np.where(t==0)[0]
        if len(t) != 0:
            t = max(t)
            assert t > 3*len(data)/4
            data[chr] = data[chr][:t]
        data[chr] = data[chr] / np.sum(data[chr])
    return data

def fit_linear_regression(data, starting_npoints = 5,
                          step = 1,
                          crop_min = -1.75,
                          crop_max = -0.05):

    maxdist = min([len(i) for i in data.values()])
    distances = np.log(np.arange(maxdist)*resolution+1)
    coeffs = []
    regression = LinearRegression()

    assert maxdist-starting_npoints-1 > 1
    starts = range(1,maxdist-starting_npoints-1, step)

    plot_distances = []
    local_average = []
    for st in starts:
        npoints = min(starting_npoints + (st*resolution)//50000, 5000000 // resolution)
        end = st + npoints
        if end >= len(distances):
            break
        X = distances[st:end].reshape(end-st,1)
        curr_coefs = []
        for chr in data:
            Y = np.log(data[chr][st:end].reshape(end-st,1))
            reg = regression.fit(X,Y)
            curr_coefs.append(reg.coef_[0][0])
        #curr_coef = np.median(curr_coefs)
        curr_coef = np.average(curr_coefs)

        # check that curr coef is not too different from last coeff
        if curr_coef < crop_min or curr_coef > crop_max:
            continue

        if (st*resolution >= 50000000):
            break

        if (st*resolution > 1000000 and len(coeffs) > 5) and False:
            # check diffference
            av = np.average(coeffs[-4:-1])
            threashold = abs(av)/2
            if abs(curr_coef - av) > threashold or len(local_average)>0: # diff is too high, probably outlayer
                local_average.append(curr_coef)
                if np.median(local_average) > 0 or np.median(local_average) < -2:
                    break
                elif abs(np.median(local_average) - av) <= threashold:
                    coeffs.append(np.median(local_average))
                    plot_distances.append(distances[st])
                    local_average = []
                continue

        coeffs.append(curr_coef)
        plot_distances.append(distances[st])

    return np.exp(plot_distances), np.array(coeffs)

def fit_delta(data, npoints = 20, step = 2):
    maxdist = min([len(i) for i in data.values()])
    distances = np.log(np.arange(maxdist)*resolution+1)
    results = []

    assert maxdist-npoints-1 > 1
    starts = range(1,maxdist-npoints*2, step)

    plot_distances = []
    for st in starts:
        curr_coefs=[]
        for chr in data:
            delta = np.average([(np.log(data[chr][i])-np.log(data[chr][i+npoints]))/np.log(data[chr][i]) for i in range(st,st+npoints)])
            curr_coefs.append(delta)
        curr_coef = np.average(curr_coefs)
        if curr_coef < 0:
            results.append(abs(curr_coef))
            plot_distances.append(distances[st])

    return np.exp(plot_distances), np.array(results)

def plot_ps(data):
    maxdist = 20000000 // resolution
    expected = data[list(data)[0]]
    for e in data.values():
        if len(e) > len(expected):
            expected = e
    distances = np.arange(1,maxdist) * resolution + 1
    return distances, expected[1:len(distances)+1]

def fit_power_low (data, npoints = 20, step = 2):
    def power_low(x, a, b):
        return  (x ** a) * b

    expected = data[list(data)[0]]
    for e in data.values():
        if len(e) > len(expected):
            expected = e

    # expected = expected / sum(expected)
    distances = np.arange(len(expected))*resolution
    coeffs = []
    # TODO think about this "20" constant
    starts = range(4,len(distances)-10,step)
    plot_distances = []
    for st in starts:
        end = st + npoints
        p0 = [-1., expected[st] * distances[st]]
        try:
            popt, pcov = curve_fit(power_low,
                               distances[st:end],
                               expected[st:end],
                               p0=p0,
                               maxfev = 1000)
        except RuntimeError: #scippy couldn't fit curve
            continue
        if 0 > popt[0] > -2:
            coeffs.append(popt[0])
            plot_distances.append(distances[st])

        # uncomment to draw fit of the curve
        # predicted = power_low(distances[st:end], *popt)
        # plt.loglog(distances[st:end], predicted)
    return plot_distances, np.array(coeffs)

def plot(X,Y,**kwargs):
    plt.semilogx(X,Y, **kwargs)

def row2color(row):
    colors={"Anopheles":"springgreen",
            "Drosophila":"mediumseagreen",
            "culex": "limegreen",
            "Aedes": "green",
            "Polypedium" : "black",
            "chick":"black",
            "mammals":"blue"
    }

    special_styles = {
            "CME" : {"linestyle":":"},
            "CIE": {"marker":"*","linestyle":":"},

            "RaoCondensinDegron": {"linestyle": "--"},
            "2017Haarhuis_KO_WAPL1" : {"linestyle":"--", "marker":  "*"},
            "2017Haarhuis_Kontrol": {"marker":  "*"},
#            "DekkerCapH-": {"color":"red","marker" : "*", "linestyle":"--"},
#            "DekkerCAPHControl": {"color": "red", "marker" : "*"},
#            "DekkerCapH2-": {"color":"salmon", "marker":"^", "linestyle" : "--"},
#            "DekkerCAPH2Control": {"color": "salmon", "marker":"^"},
#            "DekkerSMC2-":{"color": "red"},
#            "DekkerSMC2Control": {"color": "red", "linestyle":"--"},
            "DekkerPrometo":{"color":"yellow"},

            "DmelCAP":{"linestyle":"--"},
            "DmelRAD": {"linestyle": ":"},
        }

    result = {"color":colors[row.subtaxon]}
    if row["name"] in special_styles:
        for k,v in special_styles[row["name"]].items():
            result[k] = v

    return result

dataset = "datasets.csv"
juicer_tools = "juicer/juicer_tools_1.19.02.jar"
resolution = 10000
report = 170

datasets = pd.read_csv(dataset, sep="\t",
                       comment="#")
# define subplots
subplots = dict([(val,ind) for ind,val in enumerate(np.unique(datasets.taxon.values))])
print(subplots)
for ind in range(len(datasets)):
    row = datasets.iloc[ind]
    logging.info(row["name"])
    hic = row.link
    logging.info("Starting dump")
    data = dump(hic,juicer_tools,resolution)
    data = process(data)
    logging.info("Fitting...")
    # X,Y = fit_power_low(data)
    if row.taxon == "vertebrate":
        crop_min = -2.2
        crop_max = 0.1
    else:
        crop_min = -1.75
        crop_max = -0.25

    X,Y = fit_linear_regression(data, crop_min=crop_min, crop_max=crop_max)

    # X,Y = plot_ps(data)
    # X,Y = fit_delta(data)
    # X, Y = plot_ps(data)
    logging.info("Plotting...")
    plt.subplot(len(subplots),1,subplots[row.taxon]+1)
    plot(X,Y,**row2color(row))
    logging.info("Done!")
    if ind >= report and ind % report == 0:
        #plt.show()
        break
#plt.axhline(y=-1, ls="--")
plt.savefig("result.png",dpi=500)