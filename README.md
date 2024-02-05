# rerun_fwdpy => simulation_1000G_Dz_revision

## Environment setup for fwdpy11 and moments
You need to perform the following on computecanada to install fwdpy11 for simulation<br>
```
module load StdEnv/2020
module load python/3.9.6
####Minimum python>=3.9
python3 -m venv virtual_env 
source virtual_env/bin/activate
module load gsl/2.6
module load rust/1.70.0
export PATH=~/.cargo/bin:$PATH 
####This is to get cbindgen in the path 
pip install --upgrade pip 
pip install -r requirements.txt
```
If you don't want to install fwdpy11, just performing the last step is enough.<br> 

## Scripts that have been recently revised
src/simulation_moments.py <br>
src/simulation_sweep.py <br>
parameters.yaml <br>
requirements.txt <br>

## Scripts that are currently under revision
