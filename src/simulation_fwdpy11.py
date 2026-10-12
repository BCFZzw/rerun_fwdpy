import fwdpy11 
import numpy as np 
import fwdpy11.conditional_models
import msprime 
import copy


class IncompleteSweep(object):
    def __init__(self, stop_frequency: float):
        if not 0 < stop_frequency < 1:
            raise ValueError(f"stop_frequency must be in (0, 1), got {stop_frequency}")
        self.stop_frequency = stop_frequency

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
        if pop.mcounts[index] / 2 / pop.N >= self.stop_frequency:
            return fwdpy11.conditional_models.SimulationStatus.Success
        # make sure there's a valid return value
        return fwdpy11.conditional_models.SimulationStatus.Continue

def sweep_fixation(pop, params, sweep_site, seed):
    output = fwdpy11.conditional_models.selective_sweep(
            fwdpy11.GSLrng(seed),
            pop,
            params,
            sweep_site,
            fwdpy11.conditional_models.GlobalFixation(),
            return_when_stopping_condition_met = True ### stopping when sweep fixed
    )
    assert output.pop.generation == output.pop.fixation_times[0]
    print(f"Evolved {output.pop.fixation_times[0]} generation till fixation")
    return output.pop

def sweep_incomplete(pop, params, sweep_site, seed, stop_frequency):
    output = fwdpy11.conditional_models.selective_sweep(
            fwdpy11.GSLrng(seed),
            pop,
            params,
            sweep_site,
            IncompleteSweep(stop_frequency),
            return_when_stopping_condition_met = True ### stopping when sweep fixed
    )
    print(f"Evolved {output.pop.generation} generation till reaching frequency {stop_frequency}")
    return output.pop

def sweep_post_fix(pop, params, sweep_site, seed1, seed2, params_post_fix):
    output = sweep_fixation(pop, params, sweep_site, seed1)
    fixation_time = output.fixation_times[0]
    pop_post_fixed = copy.deepcopy(output)
    fwdpy11.evolvets(fwdpy11.GSLrng(seed2), pop_post_fixed, params_post_fix, simplification_interval = 100)
    print(f"Continue evolved post fixation for {pop_post_fixed.generation - fixation_time} generation")
    return pop_post_fixed

def add_neu_mutations(pop, mut_rate, seed):
    pop_with_mut = copy.deepcopy(pop)
    nmuts = fwdpy11.infinite_sites(fwdpy11.GSLrng(seed), pop_with_mut, mut_rate)
    return nmuts, pop_with_mut

def save_tree_with_mutations(pop, outfile, mut_rate, mut_seed):
    nmuts, pop_with_mut = add_neu_mutations(pop, mut_rate, mut_seed)
    print(f"{nmuts} mutations added to simulations")
    sweep_ts = pop_with_mut.dump_tables_to_tskit()
    sweep_ts.dump(outfile)
