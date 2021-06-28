###################################################
Codes for numerical analyses and simulations in
Sakamoto and Innan (2021)
###################################################

1. Folder "numerical"
Python scripts for numerical analysis are provided.
For calculations of gene loss rates conditional on L1 and L2, please use files "no_degeneration_XX.py".
It outputs the list of [R1, R2, c]

For long-term degeneration process, please use files "degeneration_XX.py".
In dominant case, it outputs the list of [L1(t), L2(t), t].
In additive case, it produces the list of [L1(t), L2(t), H(t), t].

2. Folder "simulation"
Python scripts for stochastic simulations are provided.

2-1. Folder "simulation/no_degeneration_ad_hoc"
Scripts of ad-hoc simulation are provided.

2-2. Folder "simulation/degeneration"
Scripts of simulation to investigate long-term degeneration process are provided.

2-3. Folder "simulation/variable_selection"
Scripts of simulation to investigate long-term degeneration process when fitness effect varies over genes are provided.

2-4 Folder "simulation/exact_method"
Scripts of exact simulation that is used for verifying ad-hoc method are provided.


