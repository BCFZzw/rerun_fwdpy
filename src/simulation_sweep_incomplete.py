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

class IncompleteSweep(object):
    def __call__(
        self, pop: fwdpy11.DiploidPopulation, index: int, key: tuple
    ) -> fwdpy11.conditional_models.SimulationStatus:
        if pop.mutations[index].key != key:
            # it is fixed or lost, neither of 
            # which we want
            return fwdpy11.conditional_models.SimulationStatus.Restart
        if pop.mcounts[index] == 0:
            return fwdpy11.conditional_models.SimulationStatus.Restart
        
        # Terminate the first time we see the 
        # variant get about a freq above the frequency
        if pop.mcounts[index] / 2 / pop.N >= stop_frequency:
            return fwdpy11.conditional_models.SimulationStatus.Success
        # make sure there's a valid return value
        return fwdpy11.conditional_models.SimulationStatus.Continue


def neutral_simulation(pop, params, seed):
    ### Mutations are added to the population after neutral simulations.
    pop_evolved = copy.deepcopy(pop)
    fwdpy11.evolvets(fwdpy11.GSLrng(seed), pop_evolved, params, simplification_interval = 10)
    return pop_evolved

def add_neu_mutations(pop, mut_rate, seed):
    pop_with_mut = copy.deepcopy(pop)
    nmuts = fwdpy11.infinite_sites(fwdpy11.GSLrng(seed), pop_with_mut, mut_rate)
    return nmuts, pop_with_mut

savePath = "/home/alouette/projects/ctb-sgravel/alouette/Simulation/2.Incomplete_0.95_5e5_Q10"

#Nielsen_R = 500 #4NLp
#Nielsen_theta = 0.002 #4Nmu
#Nielsen_alpha = 500 #2Ns

#rec_rate = Nielsen_R/4/n_sample/sim_region # mean # of breakpoints per diploid per generation
#mut_rate = Nielsen_theta/4/n_sample/2*sim_region #per haploid genome specified by fwdpy11
#sel_coeff = Nielsen_alpha/2/n_sample

### using scaling factors
Ne = 20000
n_sample = 2000
sim_region = int(5e5)
scaling_factor = Ne/n_sample ### scaling factor = 10
rec_rate = 1.25e-8 * scaling_factor
mut_rate = 1.44e-8 * scaling_factor
sel_coeff = 0.01 * scaling_factor
sim_gen = 200
stop_frequency = 0.95 ### global

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
            IncompleteSweep(),
            return_when_stopping_condition_met = True ### stopping when sweep fixed
        )
    except RuntimeError:
        start_seed = start_seed + 1
        continue
    print(output.pop.generation)

    ### At fixation, deep copy the tree to add mutiatons, save tree.
    nmuts, sweep_fix_mut = add_neu_mutations(output.pop, mut_rate * sim_region, seed)
    print(f"{nmuts} mutations added to sweep at fixation")
    outfile = os.path.join(savePath, "_".join(["sweep", str(seed), "incomplete", str(stop_frequency)]) + ".trees")
    sweep_ts = sweep_fix_mut.dump_tables_to_tskit()
    sweep_ts.dump(outfile)
    
    sim_i = sim_i + 1
    







