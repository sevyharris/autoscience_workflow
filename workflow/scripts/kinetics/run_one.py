# script to run with infinite time on west partition and go through one kinetics calculation
import sys
import kineticfun


reaction_index = int(sys.argv[1])
print(kineticfun.reaction_index2smiles(reaction_index))

if kineticfun.arkane_complete(reaction_index):
    print(f'Kinetics already calculated for reaction {reaction_index}')
    exit(0)

# kineticfun.run_TS_shell_calc(reaction_index, use_reverse=True)
# kineticfun.run_TS_overall_calc(reaction_index, use_reverse=True)
combos = 300

kineticfun.run_TS_shell_calc(reaction_index, max_combos=combos, max_conformers=12)
kineticfun.run_TS_overall_calc(reaction_index, max_combos=combos, max_conformers=12)
kineticfun.run_arkane_job(reaction_index)
