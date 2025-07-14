import fwdpy11 
import numpy as np 
import fwdpy11.tskit_tools
import fwdpy11.conditional_models
import sys
import os
import msprime 
import copy

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

savePath = "/home/alouette/projects/ctb-sgravel/alouette/Dz_sweep_final/simulation_review/src/test/test_simulations_post_fix_1000"

#Nielsen_R = 500 #4NLp
#Nielsen_theta = 0.002 #4Nmu
#Nielsen_alpha = 500 #2Ns

#rec_rate = Nielsen_R/4/n_sample/sim_region # mean # of breakpoints per diploid per generation
#mut_rate = Nielsen_theta/4/n_sample/2*sim_region #per haploid genome specified by fwdpy11
#sel_coeff = Nielsen_alpha/2/n_sample

### using scaling factors
Ne = 20000
n_sample = 2000
sim_region = int(1e6)
scaling_factor = Ne/n_sample ### scaling factor = 10
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
        "demography": fwdpy11.ForwardDemesGraph.tubes([n_sample], burnin=10)
}
params = fwdpy11.ModelParams(**pdict)

sweep_site = fwdpy11.conditional_models.NewMutationParameters(
    frequency=fwdpy11.conditional_models.AlleleCount(1),
    data=fwdpy11.NewMutationData(effect_size = sel_coeff, dominance=1),
    position=fwdpy11.conditional_models.PositionRange(left=(sim_region/2), right=(sim_region/2 + 1)),
)

### Post fixation parameters
### Log scale post fixation times , in terms of Ne
post_fix_time = [0.025, 0.05, 0.1, 0.25, 0.5]
post_fix_gen_arr = [int(t * n_sample) for t in post_fix_time]

simulation_iteration = 1000
start_seed = 5553
sim_i = 0


while sim_i < simulation_iteration:
    seed = start_seed + sim_i
    initial_ts = msprime.sim_ancestry(
        samples=n_sample,
        population_size=n_sample,
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

    ### At fixation, deep copy the tree to add mutiatons, save tree.
    nmuts, sweep_fix_mut = add_neu_mutations(output.pop, mut_rate * sim_region, seed)
    print(f"{nmuts} mutations added to sweep at fixation")
    outfile = os.path.join(savePath, "_".join(["sweep", str(seed), "fixation"]) + ".trees")
    sweep_ts = sweep_fix_mut.dump_tables_to_tskit()
    sweep_ts.dump(outfile)
    
    for post_fix_gen in post_fix_gen_arr:
        pdict_post_fix = {
            "recregions": [fwdpy11.PoissonInterval(0, sim_region, sim_region*rec_rate, discrete=True)],
            "gvalue": fwdpy11.Multiplicative(2.0),
            "rates": (0, 0, None),
            "prune_selected": False,
            "simlen": post_fix_gen,
            "demography": fwdpy11.ForwardDemesGraph.tubes([n_sample], burnin=10)
        }
        params_post_fix = fwdpy11.ModelParams(**pdict_post_fix)
        ### Continue neutrally evolving the population, stop at each iteration
        pop_post_fix = neutral_simulation(output.pop, params_post_fix, seed)
        assert pop_post_fix.generation == output.pop.fixation_times[0] + post_fix_gen
        ### After stopping, deep copy the popualtion and add the mutations
        nmuts, pop_with_mut = add_neu_mutations(pop_post_fix, mut_rate * sim_region, seed)
        print(f"{nmuts} mutations added to sweep at post fixation gen {post_fix_gen}")
        outfile = os.path.join(savePath, "_".join(["sweep", str(seed), "post_fix", str(post_fix_gen), "gen"]) + ".trees")
        sweep_ts = pop_with_mut.dump_tables_to_tskit()
        sweep_ts.dump(outfile)
    
    sim_i = sim_i + 1
    







