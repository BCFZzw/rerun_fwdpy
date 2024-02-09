# rerun_fwdpy => simulation_1000G_Dz_revision

## Environment setup for fwdpy11 and moments
fwdpy11 manual: https://molpopgen.github.io/fwdpy11/intro.html# <br>
moments manual: https://momentsld.github.io/moments/

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
If you only want to install moments, just remove fwdpy11 from the requirement.txt and perform the last step is enough.<br> 

## Scripts that have been recently revised
src/simulation_moments.py <br>
src/simulation_sweep.py <br>
parameters.yaml <br>
requirements.txt <br>

## Scripts that are currently under revision
pastCode/zarrMoment1000G.py <br>
This is the functions I used to compute Dz from 1000G VCFs. It's probably not in a understandable format in its current stage.<br>

### Note: 
zarrMoment1000G.py requires an input bed file that contains the SNP position of the 1000Genomes data over a pre-defined range. <br>
The pre-defined range can separate 1) different functional regions 2) job parallelization. <br>
From the region, the genotype VCF is read by scikit-allel dask function. More on scikit-allel refer to this tutorial: http://alimanfoo.github.io/2017/06/14/read-vcf.html <br>
Scikit-allel has stopped its maintenance. An alternative would be sgkit: https://pystatgen.github.io/sgkit/latest/ <br>
<br>
Prior to Moments calculation, there are 3 levels of window-parsing in various functions in zarrMoment1000G.py: <br>
1. pre-defined region <br>
2. parse the pre-defined region into 0.04 cM bins <br>
3. if the 0.04 cM bin contains too many SNPs, the memory consumption would be huge. Further subset the 0.04 cM bin into smaller bins (less than 5000 SNPs) and combine the results together. <br>

## TODO 
- [ ] The functions in zarrMoment1000G.py are composite and mashed up together. Separate into smaller and more resuable functional units. <br>
- [ ] Combine all window-parsing in a simpler chunked array function, and perform sanity check. <br>
- [ ] Follow "_" naming style instead of Camel, give clear variable and function name. Put more documentation. <br>
- [ ] Add type hinting in the functions. <br>
