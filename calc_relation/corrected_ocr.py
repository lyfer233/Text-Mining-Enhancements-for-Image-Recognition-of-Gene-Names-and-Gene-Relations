import os
from difflib import SequenceMatcher

import pandas as pd


def corrected_processing_by_dict(dictionary, ocr, thresh=0.7):
    """
    Use exHUGO dictionary to correct an ocr result by calculating similarity,
    thresh=0.7, the optimal threshold with recall,
    thresh=0.9, the optimal threshold with the mean of precision and recall,
    thresh=1.0, the optimal threshold with precision.

    :param dictionary: exHUGO dictionary, list
    :param ocr: one ocr results, string
    :param thresh: optimal threshold, ranges from 0 to 1, float
    :return:
    """
    seq_match_ratio = [SequenceMatcher(None, ocr.upper(), gene.upper()).ratio() for gene in dictionary]
    corrected_ocr = dictionary[seq_match_ratio.index(max(seq_match_ratio))] if round(max(seq_match_ratio),
                                                                                     3) >= thresh else '-'
    return corrected_ocr


if __name__ == "__main__":
    ground_truth_df = pd.read_csv("data/finalized_genes.csv")
    ground_truth_list = ground_truth_df['annotated_gene_name'].to_list()
    for root, dirs, files in os.walk(r"data/ocr_relations"):
        for file in files:
            startor_list = []
            receptor_list = []
            df = pd.read_csv(root + '/' + file)
            for index, row in df.iterrows():
                first_gene = corrected_processing_by_dict(ground_truth_list, row['startor'])
                second_gene = corrected_processing_by_dict(ground_truth_list, row['receptor'])
                if first_gene == '-' or second_gene == '-':
                    continue
                startor_list.append(first_gene)
                receptor_list.append(second_gene)

            result = {
                "startor": startor_list,
                "receptor": receptor_list,
                # "relation": df['relation'].to_list,
            }
            pd.DataFrame(result).to_csv("data/new_ocr_relations/" + file, index=False)


