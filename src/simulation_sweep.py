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
    fwdpy11.evolvets(fwdpy11.GSLrng(seed), pop, params, simplification_interval = 10)

def add_neu_mutations(pop, mut_rate, seed):
    pop_with_mut = copy.deepcopy(pop)
    nmuts = fwdpy11.infinite_sites(fwdpy11.GSLrng(seed), pop_with_mut, mut_rate)
    return nmuts, pop_with_mut

savePath = "/home/alouette/projects/ctb-sgravel/alouette/Dz_sweep_final/simulation_review/src/test/test_simulations"

Nielsen_R = 500 #4NLp
Nielsen_theta = 0.002 #4Nmu
Nielsen_alpha = 500 #2Ns

n_sample = 500
sim_region = int(1e6)
rec_rate = Nielsen_R/4/n_sample/sim_region # mean # of breakpoints per diploid per generation
mut_rate = Nielsen_theta/4/n_sample/2*sim_region #per haploid genome specified by fwdpy11
sel_coeff = Nielsen_alpha/2/n_sample

pdict = {
        "recregions": [fwdpy11.PoissonInterval(0, sim_region, sim_region*rec_rate, discrete=True)],
        # Here, gvalue as multiplicative(2.0) means 1, 1+hs, 1+2s.
        "gvalue": fwdpy11.Multiplicative(2.0),
        "rates": (0, 0, None),
        "prune_selected": False,
        "simlen": 200,
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
post_fix_gen = 10
post_fix_iteration = 10
pdict_post_fix = pdict.copy()
pdict_post_fix["simlen"] = post_fix_gen
params_post_fix = fwdpy11.ModelParams(**pdict_post_fix)

simulation_iteration = 100
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
    ### Record fixation time

    for i in range(0, post_fix_iteration):
        neutral_simulation(output.pop, params_post_fix, seed)
        nmuts, pop_with_mut = add_neu_mutations(output.pop, mut_rate, seed)
        print(f"{nmuts} mutations added to sweep at post fixation gen {i*post_fix_gen}")
        outfile = os.path.join(savePath, "_".join(["sweep", str(seed), "post_fix", str(i*post_fix_gen), "gen"]) + ".trees")
        sweep_ts = pop_with_mut.dump_tables_to_tskit()
        sweep_ts.dump(outfile)
    
    sim_i = sim_i + 1

### what time I put for neutral evolution?
#neutral_pop = fwdpy11.DiploidPopulation(n_sample, sim_region)
#neutral_simulation(neutral_pop, params, seed)
#nmuts, pop_with_mut = add_neu_mutations(neutral_pop, mut_rate, seed)
#print(f"{nmuts} mutations added to neutral background")
#neutral_prefix = os.path.join(savePath, "neutral", "neutral_" + str(seed))
#neutral_ts = neutral_pop.dump_tables_to_tskit()






