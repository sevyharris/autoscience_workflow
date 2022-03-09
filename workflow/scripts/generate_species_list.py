import pandas as pd
import cantera as ct
import rmgpy.chemkin
import rmgpy.species


# model_cti = '/work/westgroup/harris.se/autoscience/autoscience_workflow/resources/nheptane_1.cti'
chemkin_path = "/work/westgroup/harris.se/autoscience/autoscience_workflow/resources/nheptane1/chem.inp"
dictionary_path = "/work/westgroup/harris.se/autoscience/autoscience_workflow/resources/nheptane1/species_dictionary.txt"
transport_path = "/work/westgroup/harris.se/autoscience/autoscience_workflow/resources/nheptane1/tran.dat"
species_list, reaction_list = rmgpy.chemkin.load_chemkin_file(
    chemkin_path,
    dictionary_path=dictionary_path,
    transport_path=transport_path
)

print(species_list[0].molecule[0].smiles)

entries = []
for i, species in enumerate(species_list):
    entries.append([i, str(species), species.molecule[0].smiles])

df = pd.DataFrame(entries, columns=['i', 'name', 'SMILES'])
df.to_csv('species_list.csv')

#gas = ct.Solution(model_cti)

#print(gas.species_names)
#species_list = gas.species()

#print(species_list[3])

