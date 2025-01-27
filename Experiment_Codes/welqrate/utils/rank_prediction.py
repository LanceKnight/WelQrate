# Open the input file for reading

def rank_prediction(type, model_name):
    with open(f'result/{model_name}/per_molecule_pred_of_{type}_set.txt', 'r') as f:
        data = [(float(line.split('\t')[0]), line.split('\t')[1] ) for line in
                f.readlines()]

    ranked_data = sorted(data, key=lambda x: x[0], reverse=True)

    with open(f'result/{model_name}/ranked_mol_score_{type}.txt', 'w') as f:
        for i, (score, label) in enumerate(ranked_data):
            f.write(f"{i}\t{score}\t{label}")