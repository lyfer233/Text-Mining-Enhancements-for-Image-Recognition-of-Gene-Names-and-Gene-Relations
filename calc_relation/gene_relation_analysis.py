import os
import pandas as pd


def calc(will_check_list, ground_truth_list):
    precision = 0.0
    recall = 0.0
    overleap_count = 0
    for gene_pair1 in will_check_list:
        for gene_pair2 in ground_truth_list:
            if (gene_pair1[0] == gene_pair2[0] and gene_pair1[1] == gene_pair2[1]) or (
                    gene_pair1[0] == gene_pair2[1] and gene_pair1[1] == gene_pair2[0]):
                overleap_count += 1
                break

    # If will_check_list is empty, then it cannot be divided by 33，
    if (len(will_check_list) != 0):
        precision = 1.0 * overleap_count / len(will_check_list)
    if (len(ground_truth_list) != 0):
        recall = 1.0 * overleap_count / len(ground_truth_list)

    return precision, recall


if __name__ == '__main__':
    # original data
    filename_pmcid_pair_file = 'data/figurename.xlsx'
    ground_truth_gene_relations_file = '45_gene_relation_annotation.csv'
    ocr_results_folder = 'data/new_ocr_relations/'
    co_occur_same_sent_folder = 'data/co_occur_same_sent/'
    co_occur_neigh_folder = 'data/co_occur_neigh_sent/'
    openie_folder = 'data/openie/'
    full_text_path = "../hefei_pipeline/calc_indicator/GC_with_meta_results_GT_yijie/co_occur_caption_sent/"

    # figure name & corresponding pmcid xlsx
    df_filename_pmcid_pair = pd.read_excel(filename_pmcid_pair_file)
    cols = ['filename', 'PMCID']
    df_filename_pmcid_pair[cols] = df_filename_pmcid_pair[cols].apply(lambda x: x.str.strip())
    figure_names = list(df_filename_pmcid_pair['filename'])

    # ground truth csv -- gene relations
    df_gene_relations_ground_truth = pd.read_csv(ground_truth_gene_relations_file)

    fig_name = []
    ocr_in_total = []
    ground_truth_in_total = []
    ocr_vs_ground_truth_overlap = []
    same_sent_ocr_vs_ground_truth_overlap = []
    neigh_sent_ocr_vs_pubtator_overlap = []
    openie_ocr_vs_ground_truth_overlap = []
    my_result = []

    total_precision = 0.0
    total_recall = 0.0
    for jpg_file in figure_names:
        fig_name.append(jpg_file)
        jpg_file_name = os.path.basename(jpg_file[:-4])
        pmcid = jpg_file_name.split('_')[0]
        ocr_relations_file_name = ocr_results_folder + pmcid + '_ocr_results.csv'
        same_sent_file_name = co_occur_same_sent_folder + pmcid + '_co_occurrence_same_sent.csv'
        neigh_sent_file_name = co_occur_neigh_folder + pmcid + '_co_occurrence_neighboring.csv'
        openie_file_name = openie_folder + pmcid + '_openie.csv'
        my_result_file_name = full_text_path + pmcid + ".csv"

        # ground truth gene relations, gene names in each gene relation
        df_ground_truth_gene_relations = df_gene_relations_ground_truth.loc[
            df_gene_relations_ground_truth['fig_name'] == jpg_file]
        ground_truth_in_total.append(str(len(df_ground_truth_gene_relations.index)))
        ground_truth_tuple_list = []
        for index, row in df_ground_truth_gene_relations.iterrows():
            ground_truth_tuple_list.append((row['activator'], row['receptor']))

        print("{},{}".format(pmcid, len(ground_truth_tuple_list)))

        # ocr relations
        df_ocr = pd.read_csv(ocr_relations_file_name)
        ocr_in_total.append(str(len(df_ocr.index)))
        ocr_tuple_list = []
        for index, row in df_ocr.iterrows():
            ocr_tuple_list.append((row['startor'], row['receptor']))

        precision, recall = calc(ocr_tuple_list, ground_truth_tuple_list)
        total_precision += precision
        total_recall += recall

        # gene co-occurence relations -- same
        df_same = pd.read_csv(same_sent_file_name)
        same_sent_tuple_list = []
        # same_sent_tuple_list += ocr_tuple_list
        for index, row in df_same.iterrows():
            same_sent_tuple_list.append((row['gene_name_1'], row['gene_name_2']))
        # precision, recall = calc(same_sent_tuple_list, ground_truth_tuple_list)
        # total_precision += precision
        # total_recall += recall

        # gene co-occurence relations -- neighbor
        df_neigh = pd.read_csv(neigh_sent_file_name)
        neigh_sent_tuple_list = []
        # neigh_sent_tuple_list += ocr_tuple_list
        for index, row in df_neigh.iterrows():
            neigh_sent_tuple_list.append((row['gene_name_1'], row['gene_name_2']))
        # precision, recall = calc(neigh_sent_tuple_list, ground_truth_tuple_list)
        # total_precision += precision
        # total_recall += recall

        # openie extracted relations
        df_openie = pd.read_csv(openie_file_name)
        openie_tuple_list = []
        openie_tuple_list += ocr_tuple_list
        for index, row in df_openie.iterrows():
            if row['confidence_score'] == 1:
                openie_tuple_list.append((row['subject'], row['object']))

        # # my result
        # try:
        #     full_text_df = pd.read_csv(my_result_file_name)
        # except Exception as e:
        #     my_result.append(0)
        #     continue
        # full_text_list = []
        # full_text_list += ocr_tuple_list
        # for index, row in full_text_df.iterrows():
        #     full_text_list.append((str(row['gene_name_1']).upper(), str(row['gene_name_2']).upper()))

        ocr_count = 0
        same_sent_count = 0
        neigh_sent_count = 0
        openie_count = 0
        my_result_count = 0

        for tuple_ground_truth in ground_truth_tuple_list:

            # # my result vs ground truth
            # my_flag = 0
            # for full_text_tuple in full_text_list:
            #     if my_flag == 0 and full_text_tuple[0] == tuple_ground_truth[0] and full_text_tuple[1] == tuple_ground_truth[1]:
            #         my_result_count += 1
            #         my_flag = 1
            #     elif my_flag == 0 and full_text_tuple[0] == tuple_ground_truth[1] and full_text_tuple[1] == tuple_ground_truth[0]:
            #         my_result_count += 1
            #         my_flag = 1
            # ocr vs ground truth
            ocr_flag = 0
            for tuple_ocr in ocr_tuple_list:
                if ocr_flag == 0 and tuple_ocr[0] == tuple_ground_truth[0] and tuple_ocr[1] == tuple_ground_truth[1]:
                    ocr_count += 1
                    ocr_flag = 1
                elif ocr_flag == 0 and tuple_ocr[1] == tuple_ground_truth[0] and tuple_ocr[0] == tuple_ground_truth[1]:
                    ocr_count += 1
                    ocr_flag = 1

            # same_sent vs ground truth
            same_sent_flag = 0
            for tuple_same_sent in same_sent_tuple_list:
                if same_sent_flag == 0 and tuple_same_sent[0] == tuple_ground_truth[0] and tuple_same_sent[1] == \
                        tuple_ground_truth[1]:
                    same_sent_count += 1
                    same_sent_flag = 1
                elif same_sent_flag == 0 and tuple_same_sent[1] == tuple_ground_truth[0] and tuple_same_sent[0] == \
                        tuple_ground_truth[1]:
                    same_sent_count += 1
                    same_sent_flag = 1

            # neigh_sent vs ground truth
            neigh_sent_flag = 0
            for tuple_neigh_sent in neigh_sent_tuple_list:
                if neigh_sent_flag == 0 and tuple_neigh_sent[0] == tuple_ground_truth[0] and tuple_neigh_sent[1] == \
                        tuple_ground_truth[1]:
                    neigh_sent_count += 1
                    neigh_sent_flag = 1
                elif neigh_sent_flag == 0 and tuple_neigh_sent[1] == tuple_ground_truth[0] and tuple_neigh_sent[0] == \
                        tuple_ground_truth[1]:
                    neigh_sent_count += 1
                    neigh_sent_flag = 1

            # openie vs ground truth
            openie_flag = 0
            for tuple_openie in openie_tuple_list:
                if openie_flag == 0 and tuple_ground_truth[0] in tuple_openie[0] and tuple_ground_truth[1] in \
                        tuple_openie[1]:
                    openie_count += 1
                    openie_flag = 1
                elif openie_flag == 0 and tuple_ground_truth[0] in tuple_openie[1] and tuple_ground_truth[1] in \
                        tuple_openie[0]:
                    openie_count += 1
                    openie_flag = 1

        ocr_vs_ground_truth_overlap.append(str(ocr_count))
        same_sent_ocr_vs_ground_truth_overlap.append(str(same_sent_count))
        neigh_sent_ocr_vs_pubtator_overlap.append(str(neigh_sent_count))
        openie_ocr_vs_ground_truth_overlap.append(str(openie_count))
        # my_result.append(str(my_result_count))

    dict = {
        'fig_name': fig_name,
        'ocr_in_total': ocr_in_total,
        'ground_truth_in_total': ground_truth_in_total,
        'ocr_vs_ground_truth_overlap': ocr_vs_ground_truth_overlap,
        'same_sent_ocr_vs_ground_truth_overlap': same_sent_ocr_vs_ground_truth_overlap,
        'neigh_sent_ocr_vs_pubtator_overlap': neigh_sent_ocr_vs_pubtator_overlap,
        'openie_ocr_vs_ground_truth_overlap': openie_ocr_vs_ground_truth_overlap
    }

    print(total_precision/33.0)
    print(total_recall/33.0)
    df = pd.DataFrame(dict)
    # df.to_csv('data/ocr_updated.csv', index=False)
