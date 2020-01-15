# QCEIMS-analysis
1. deal with different kinds of file related to QCEIMS and NIST 17, for example, sdf, msp, jdx,mol;

2. ways to calculate scores including Cos, Dot and Match similarity, according to:

   https://doi.org/10.1016/1044-0305(94)87009-8, Stein and Scott 1994
   https://doi.org/10.1016/1044-0305(94)85022-4, Stein 1994
3. plot functions (MSpec.py)

4. qceims.out analysis is done by qceimsout.py, which has two modes: -f and -b 

## jobarray tool:

I wrote this small tool to run array through slurm system.
It is flexible with CPUs and can help you avoid queuing!
To use this tool:

1)copy the array.bash and bin file and the tmol file to your folder. 
2)You need to change the settings in array.bash code, 
	including the 
	$workdir and $user
3)run 'bash array.bash [STUCTURE].tmol’

4)after the jobs finished, you can find all the files including result.jdx 
Under $workdir/[STUCTURE] folder.
