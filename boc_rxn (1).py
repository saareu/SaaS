from rdkit import Chem
from rdkit.Chem import AllChem

def perform_reaction(head_smiles, tail_smiles, reaction_smirks):
    # Convert the SMILES strings to RDKit molecules
    head_molecule = Chem.MolFromSmiles(head_smiles)
    tail_molecule = Chem.MolFromSmiles(tail_smiles)

    # Create the reaction from the SMIRKS string
    reaction = AllChem.ReactionFromSmarts(reaction_smirks)

    # Run the reaction
    reaction_result = reaction.RunReactants((head_molecule, tail_molecule))

    # Check if the reaction produced any product
    if reaction_result:
        product = reaction_result[0][0]  # Take the first product (assuming one product is formed)
        return Chem.MolToSmiles(product)
    else:
        return None

def sequential_reaction(head_smiles, reactions):
    # Perform the reaction sequentially for each tail
    for i, (tail_smiles, repeat) in enumerate(reactions, start=1):
        print(f"Tail {i}:")
        for _ in range(repeat):
            head_smiles = perform_reaction(head_smiles, tail_smiles, reaction_smirks)
            if not head_smiles:
                print("No products formed.")
                return
        print(f"Product after {repeat} reactions: {head_smiles}\n")

if __name__ == "__main__":
    # Get inputs from the user
    head_smiles = input("Enter the SMILES for the HEAD compound: ")
    num_tails = int(input("Enter the number of different tails: "))

    # Collect SMILES and repeat count for each tail
    reactions = []
    for i in range(num_tails):
        tail_smiles = input(f"Enter the SMILES for tail {i+1}: ")
        repeat = int(input(f"Enter the number of repeats for tail {i+1}: "))
        reactions.append((tail_smiles, repeat))

    # Get SMIRKS reaction pattern from the user
    reaction_smirks = input("Enter the SMIRKS reaction: ")

    # Perform the sequential reaction
    sequential_reaction(head_smiles, reactions)
