import fwdpy11 
import numpy as np 
import fwdpy11.conditional_models
import os
import msprime 
import copy
import argparse


### Update: directly save the tree sequences, without converting to VCFs.
### Read tree sequences using tskit for moments.

parser = argparse.ArgumentParser(
                    prog='Dz_simulation_fixation',
                   )
parser.add_argument('-s', '--savedir', dest = "savedir", required = True)
parser.add_argument('-f', '--sampling-factor', dest = "sampling_factor", help = "How much to sample from output.", required = False, default = 1)
parser.add_argument('-p', '--post-fix', dest = "post_fix_gen", help = "How much time post fixation, scaled to simulations.", required = False, default = 0)
args = parser.parse_args()


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
    sample_ind_ids = rng.choice(individual_ids, size = n_sample, replace = False)
    ### extract nodes for diploid individuals
    keep_nodes = []
    for i in sample_ind_ids:
        keep_nodes.extend(tree.individual(i).nodes)
    sampled_tree = tree.simplify(keep_nodes)
    return sampled_tree




savePath = args.savedir
sampling_factor = float(args.sampling_factor)
post_fix_gen = int(args.post_fix_gen)

os.makedirs(savePath, exist_ok=True)

#Nielsen_R = 500 #4NLp
#Nielsen_theta = 0.002 #4Nmu
#Nielsen_alpha = 500 #2Ns

#rec_rate = Nielsen_R/4/n_sample/sim_region # mean # of breakpoints per diploid per generation
#mut_rate = Nielsen_theta/4/n_sample/2*sim_region #per haploid genome specified by fwdpy11
#sel_coeff = Nielsen_alpha/2/n_sample

### using scaling factors
Ne = 20000
ne_scaled= 2000
scaling_factor = Ne/ne_scaled ### scaling factor = 10
rec_rate = 1.25e-8 * scaling_factor
mut_rate = 1.44e-8 * scaling_factor
sel_coeff = 0.01 * scaling_factor
n_sample = int(ne_scaled * sampling_factor)

sim_region = int(5e5)
sim_gen = 200

simulation_iteration = 1000
start_seed = 5553
sim_i = 0
rng = np.random.default_rng(seed=42)



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

sweep_site = fwdpy11.conditional_models.NewMutationParameters(
    frequency=fwdpy11.conditional_models.AlleleCount(1),
    data=fwdpy11.NewMutationData(effect_size = sel_coeff, dominance=1),
    position=fwdpy11.conditional_models.PositionRange(left=(sim_region/2), right=(sim_region/2 + 1)),
)


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
    ### A catch to cath failed tree spans due to random recombination breakpoints, continue with next simulations
    try:
        output = fwdpy11.conditional_models.selective_sweep(
            fwdpy11.GSLrng(seed),
            pop,
            params,
            sweep_site,
            fwdpy11.conditional_models.GlobalFixation(),
            return_when_stopping_condition_met = True ### stopping when sweep fixed
        )
    except RuntimeError:
        start_seed = start_seed + 1
        continue
    assert output.pop.generation == output.pop.fixation_times[0]

    pop = copy.deepcopy(output.pop)
    outfile = os.path.join(savePath, "_".join(["sweep", str(seed), "fixation"]) + ".trees")
    if post_fix_gen != 0:
        outfile = os.path.join(savePath, "_".join(["sweep", str(seed), "post_fix", str(post_fix_gen), "gen"]) + ".trees")
        pdict_post_fix = {
            "recregions": [fwdpy11.PoissonInterval(0, sim_region, sim_region*rec_rate, discrete=True)],
            "gvalue": fwdpy11.Multiplicative(2.0),
            "rates": (0, 0, None),
            "prune_selected": False,
            "simlen": post_fix_gen,
            "demography": fwdpy11.ForwardDemesGraph.tubes([ne_scaled], burnin=10)
        }
        params_post_fix = fwdpy11.ModelParams(**pdict_post_fix)
        ### Continue neutrally evolving the population, stop at each iteration
        pop = neutral_simulation(output.pop, params_post_fix, seed)
        assert pop.generation == output.pop.fixation_times[0] + post_fix_gen
        


    nmuts, pop_with_mut = add_neu_mutations(pop, mut_rate * sim_region, seed)
    print(f"{nmuts} mutations added to sweep at post fixation gen {post_fix_gen}")
    
    sweep_ts = pop_with_mut.dump_tables_to_tskit()
    sampled_ts = sampling_individuals(sweep_ts, n_sample)
    sampled_ts.dump(outfile)
    
    sim_i = sim_i + 1
    







