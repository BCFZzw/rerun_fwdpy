import tskit
import numpy as np
import os
import glob
from simulation_selection_detection import *
import argparse


parser = argparse.ArgumentParser(
                    prog='Simulation_LD_calculation',
                    description='place_holder'
                   )
parser.add_argument('-d', '--dataPath', dest = "dataPath", required = True)
parser.add_argument('-s', '--savePath', dest = "savePath", required = False, default = ".")
parser.add_argument('-w', '--windows', dest = "n_wins", required = True)
parser.add_argument('-L', '--length', dest = "length", help = "The length of region simulated.", required = True)
parser.add_argument('-n', '--name', dest = "name", help = "The output file name", required = True)
parser.add_argument('--suffix', dest = "suffix", help = "Specify which file to parse LD by specifying the suffix. Default is .trees.", required = False, default = ".trees")
args = parser.parse_args()


dataPath = args.dataPath
savePath = args.savePath
n_wins = int(args.n_wins)
sim_region = int(args.length)
name = args.name
suffix = args.suffix
window_size = sim_region/n_wins
windows = [(i*window_size, (i+1)*window_size) for i in range(n_wins)]

os.makedirs(savePath, exist_ok=True)

filelist = glob.glob(os.path.join(dataPath, "*" + suffix))
filelist.sort()

D2_arr = []
Dz_arr = []
pi2_arr = []
for infile in filelist:
        base = infile.split("/")[-1].split(".")[0]
        ts = tskit.load(os.path.join(dataPath, infile))
        genotype, pos_array = tskit_to_allel(ts)
        D2, Dz, pi2 = window_LD(genotype, pos_array, windows)
        D2_arr.append(D2)
        Dz_arr.append(Dz)
        pi2_arr.append(pi2)


### combine all LD data into one
LD_dict = {"D2": D2_arr, "Dz": Dz_arr, "pi2": pi2_arr}
np.save(os.path.join(savePath, "_".join(["LD_dict", str(sim_region), str(n_wins), "all", name + ".npy"])), LD_dict)