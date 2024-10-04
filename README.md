# BAI-in-PSLB
This is the code for the paper ``Almost Minimax Optimal Best Arm Identification in Piecewise Stationary Linear Bandits''

## Description of the codes
- Code for G-optimal allocation computation from [Minimum-Volume Ellipsoids: Theory and Algorithms](http://archive.siam.org/books/mo23/):
	- minvol.m
	- initwt.m
	- updateR.m
	- updatevar.m
	- ellipse.m
- PS$\varepsilon$BAI$^+$ and PS$\varepsilon$BAI algorithms and their subroutines:
	- PSBAI_plus_func.m
	- PSBAI_func.m
	- LCD.m
	- LCA.m
	- rho_func.m
	- EUplus.m
- The strong baselines:
	- DBAI.m
	- DBAIbeta.m
- Functions related to the instance:
	- sample_n_arm.m
	- ClipG.m
	- context_generator.m
	- Fixed_CP_generation.m
- To get Figure 2 and 3, run the following scripts in order (the first two scripts can be skipped as the results are provided in "trial_1.m", "trial_2.m","trial_3.m"):
	- RUN_THIS.m: conduct PS$\varepsilon$BAI$^+$ on the 12 instances 20 times
	- RUN_THIS_PSBAI.m: conduct PS$\varepsilon$BAI on the 12 instances 20 times with misspecified $L_{\max}$
	- DBAI.m
	- DBAIbeta.m
	- draw.m
- To get Figure 4, i.e., run the following scripts in order (the first two scripts can be skipped as the results are provided in "trial_2.m", "trial_gamma_1.mat.m", "trial_gamma_2.mat.m","trial_gamma_4.mat.m"):
	- RUN_THIS_PSBAI_gamma.m
	- draw_ablation_gamma.m
