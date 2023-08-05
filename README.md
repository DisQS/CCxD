# CCxD
Codes to simulate the real-space RG in variants of the Chalker-Coddington models. Two separate codes exist, CC2D reproduces the method found in [Real-Space renormalization group approach to the integer quantum Hall effect](https://arxiv.org/abs/cond-mat/0501246). CCGD has increased functionality to implement geometric disorder to the RG scheme. CCGD at 0 GD can reproduce CC2D; all jobs that could be run on CC2D could also be run on CCGD.

## RG scripts
Within each folder there are two main shell scripts, a master.sh and pushmaster.sh. Each script submits multiple jobs via a slurm queue and handles its own paralellisation via slurm dependencies. Each iteration of the code consists of paralellisable computation, followed by averaging over distributed jobs, and preparation for the next iteration.

### master.sh
The master script combines multiple functionalities from earlier versions of the code, allowing one single code run to create a data set. In this case the master script is used to find the fixed point distribution. There are a number of flags submitted with the job execution to specify the parameters of the individual job. These are as follows:

   `-c`: number of configurations, the default is 1M but normally 100M would be run  

   `-m`: matrix type used, either specify 'C' or 'S' to use a S-matrix with a basis of t or with theta, respectively  

   `-i`: the amount of iterations to run the job for.  

   `-s`: no argument is passed, this specifies whether the distribution in terms of saddle point heights is symmetrised after each iteration  

   `-b`: for better bounds on error estimation, the job can be run multiple times simultaneously, this parameter specifies how many times (usually 5)  

   `-r`: to restart the process from a previous point. The argument is the iteration number that you want to start from  

   `-n`: once restarted, this specifies the number of iterations from -r you would like to iterate for  

   `-g`: the proportion of geometric disorder introduced to the system  

With this in mind a sample job running command could look like the following  



   ```./master.sh -c 100000000 -s -i 10 -b 5 -m C -g 0.1```
which would run the code for 10 iterations, 5 times over, on 100M configurations with the matrix in the t basis, including p=0.1 of geometric disorder  

If this does not result in a FP distribution the code can be run again as follows  


`./master.sh -c 100000000 -s -i 10 -b 5 -m C -g 0.1 -r 10 -n 10`
This would run the code for another 10 iterations.