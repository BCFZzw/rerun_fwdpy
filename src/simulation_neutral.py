import fwdpy11 
import numpy as np 
import fwdpy11.tskit_tools
import fwdpy11.conditional_models
import sys
import os
import msprime 
import copy
import math

### Adapted from fwdpy11 manual on selective sweep, https://molpopgen.github.io/fwdpy11/short_vignettes/incomplete_sweep.html


### Update: directly save the tree sequences, without converting to VCFs.
### Read tree sequences using tskit for moments.



def neutral_simulation(pop, params, seed):
    ### Mutations are added to the population after neutral simulations.
    pop_evolved = copy.deepcopy(pop)
    fwdpy11.evolvets(fwdpy11.GSLrng(seed), pop_evolved, params, simplification_interval = 10)
    return pop_evolved

def add_neu_mutations(pop, mut_rate, seed):
    pop_with_mut = copy.deepcopy(pop)
    nmuts = fwdpy11.infinite_sites(fwdpy11.GSLrng(seed), pop_with_mut, mut_rate)
    return nmuts, pop_with_mut


def sampling_individuals(tree, n_sample):
    n_individuals = tree.num_individuals
    if n_sample == n_individuals:
        return copy.deepcopy(tree)
    individual_ids = list(range(n_individuals)) ### currently it's just this list in 1 pop simulation
    sample_ind_ids = np.random.choice(individual_ids, size = n_sample, replace = False)
    ### extract nodes for diploid individuals
    keep_nodes = []
    for i in sample_ind_ids:
        keep_nodes.extend(tree.individual(i).nodes)
    sampled_tree = tree.simplify(keep_nodes)
    return sampled_tree




savePath = "/home/alouette/projects/ctb-sgravel/alouette/Simulation/0.neutrality"

#Nielsen_R = 500 #4NLp
#Nielsen_theta = 0.002 #4Nmu
#Nielsen_alpha = 500 #2Ns

#rec_rate = Nielsen_R/4/n_sample/sim_region # mean # of breakpoints per diploid per generation
#mut_rate = Nielsen_theta/4/n_sample/2*sim_region #per haploid genome specified by fwdpy11
#sel_coeff = Nielsen_alpha/2/n_sample

### using scaling factors
Ne = 20000
ne_scaled= 2000
sampling_factor = 1
n_sample = int(ne_scaled * sampling_factor)
sim_region = int(5e5)
scaling_factor = Ne/ne_scaled ### scaling factor = 10
rec_rate = 1.25e-8 * scaling_factor
mut_rate = 1.44e-8 * scaling_factor
sel_coeff = 0.01 * scaling_factor
sim_gen = 200

pdict = {
        "recregions": [fwdpy11.PoissonInterval(0, sim_region, sim_region*rec_rate, discrete=True)],
        # Here, gvalue as multiplicative(2.0) means 1, 1+hs, 1+2s.
        "gvalue": fwdpy11.Multiplicative(2.0),
        "rates": (0, 0, None),
        "prune_selected": False,
        "simlen": sim_gen,
        ### burnin minimum 10 from heuristic and from Ferrari et al, 2025
        "demography": fwdpy11.ForwardDemesGraph.tubes([ne_scaled], burnin=10)
}
params = fwdpy11.ModelParams(**pdict)


simulation_iteration = 1000
start_seed = 5553
sim_i = 0


while sim_i < simulation_iteration:
    seed = start_seed + sim_i
    initial_ts = msprime.sim_ancestry(
        samples=ne_scaled,
        population_size=ne_scaled,
        recombination_rate=rec_rate,
        random_seed=seed,
        sequence_length=sim_region,
    )
    pop = fwdpy11.DiploidPopulation.create_from_tskit(initial_ts)
    
    neutral_pop = neutral_simulation(pop, params, seed)
    nmuts, neutral_pop_mut = add_neu_mutations(neutral_pop, mut_rate * sim_region, seed)
    print(f"{nmuts} mutations added to neutral simulations")
    outfile = os.path.join(savePath, "_".join(["neutral", str(seed)]) + ".trees")
    neutral_ts = neutral_pop_mut.dump_tables_to_tskit()
    sampled_ts = sampling_individuals(neutral_ts, n_sample)
    sampled_ts.dump(outfile)
    
    sim_i = sim_i + 1
    







