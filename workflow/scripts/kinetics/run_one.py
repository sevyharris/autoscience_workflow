# script to run with infinite time on west partition and go through one kinetics calculation
import sys
import kineticfun


reaction_index = int(sys.argv[1])
print(kineticfun.reaction_index2smiles(reaction_index))

kineticfun.run_TS_shell_calc(reaction_index)
kineticfun.run_TS_overall_calc(reaction_index)
