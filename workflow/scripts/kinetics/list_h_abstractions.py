import kineticfun

import autotst.species
import autotst.reaction
import autotst.calculator.gaussian


N = kineticfun.get_num_reactions()
h_abstractions = []
for reaction_index in range(0, N):
    reaction_smiles = kineticfun.reaction_index2smiles(reaction_index)
    if "[C-]#[O+]" in reaction_smiles:  # might be able to get by with replacing it as [C]#[O]
        print(f'skipping {reaction_smiles}')
        continue

    # print(reaction_smiles)
    reaction = autotst.reaction.Reaction(label=reaction_smiles)
    try:
        reaction.get_labeled_reaction()
    except AssertionError:
        print(f'skipping {reaction_smiles}')
        continue
    print(reaction_smiles, '\t', reaction.rmg_reaction.family)
    h_abstractions.append(f'{reaction_index}:' + '\t' + f'{reaction_smiles}' + '\t' + f'{reaction.rmg_reaction.family}')

for reaction in h_abstractions:
    print(reaction)
