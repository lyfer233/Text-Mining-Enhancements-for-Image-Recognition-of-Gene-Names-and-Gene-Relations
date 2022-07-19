import glob
import os
import json
import xml.etree.ElementTree as ET
import pandas as pd
from fuzzywuzzy import fuzz, process
from fuzzywuzzy.process import default_processor
import nltk
from collections import OrderedDict
from subprocess import *
import re
import itertools


# use underline '_' to replace space ' ' for each gene, to make each gene name consisting of 'one word', then tokenize
# each sentence in full text
def preprocess_sent_list_and_gene_list(full_text, gene_list):
    '''
    :param full_text: str, full text of a PMCID article
    :param gene_list: list, a list of gene names
    :return: 2 lists, processed_sent_list -> list, with each sentence as a string linked by '_'
                      processed gene_list -> list, use underline '_' to replace space ' ' for each gene
    '''

    processed_gene_list = []
    for gene in gene_list:
        gene = gene.upper()
        if ' ' in gene:
            gene_underline = gene.replace(' ', '_')
            processed_gene_list.append(gene_underline)
            full_text = full_text.replace(gene, gene_underline)
        else:
            processed_gene_list.append(gene)

    sent_list = nltk.sent_tokenize(full_text)
    processed_sent_list = []
    for sent in sent_list:
        processed_sent_list.append('_'.join(nltk.word_tokenize(sent)).upper())
    # This is an apple
    # This_is_an_apple
    # This is the (EKT) gene
    # This_is_the_(_EKT_)_gene
    # 区别1加了_
    return processed_sent_list, list(set(processed_gene_list))


def parse_json(json_file_path, pmcid, output_folder):
    '''
    :param json_file_path: str, json file contains relations extracted from figures
    :return: startor_receptor_relation_tuple_list: tuple list, contains startor-receptor-relation tuple
    '''

    with open(json_file_path) as f:
        data = json.load(f)

    startor_receptor_relation_tuple_list = []

    startor = []
    receptor = []
    relation = []
    for key, value in data.items():
        if data[key]['startor'] == None or data[key]['receptor'] == None:
            continue

        startor_receptor_relation_tuple_list.append((data[key]['startor'].strip(), data[key]['receptor'].strip(), data[key]['relation_category']))
        startor.append(data[key]['startor'].strip())
        receptor.append(data[key]['receptor'].strip())
        relation.append(data[key]['relation_category'])

    dict = {'startor': startor,
            'receptor': receptor,
            'relation': relation}

    df = pd.DataFrame(dict)
    df.to_csv(output_folder + pmcid + '_ocr_results.csv', index=False)

    return startor_receptor_relation_tuple_list


def co_occurrence_same_sentence_match(sentence_list, gene_list, output_folder):
    '''
     :param sentence_list: nested list, a list of sentences with word tokenized
     :param gene_name_set: set, a set of gene names
     :param csv_folder_path: str, a folder path to store gene co-occurrence csv file
     :param pmcid: str, pmcid, refer to the output gene co-occurrence csv file name
     :return: gene_co_occurrence_tuple_list:
                      a tuple list, which contains gene co-occurrence information
                      gene_name_1, gene_name_2, co_occurrence
     '''

    occurrences = OrderedDict((name, OrderedDict((name, 0) for name in gene_list)) for name in gene_list)

    # Find the co-occurrences:
    gene_set = set(gene_list)
    for sent in sentence_list:
        intersection = gene_set.intersection(set(sent))
        if len(intersection) > 1:
            name_list = list(intersection)
            for i in range(len(name_list)):
                for item in name_list[:i] + name_list[i + 1:]:
                    occurrences[name_list[i]][item] += 1

    # print(' ', ' '.join(occurrences.keys()))
    gene_name_1 = []
    gene_name_2 = []
    gene_co_occurrence = []
    for name_1, values in occurrences.items():
        for name_2, co_occurrence in values.items():
            if (co_occurrence != 0) and not (name_2 in gene_name_1 and name_1 in gene_name_2):
                gene_name_1.append(name_1)
                gene_name_2.append(name_2)
                gene_co_occurrence.append(co_occurrence)

    gene_name_1 = [name.replace('_', ' ').upper() for name in gene_name_1]
    gene_name_2 = [name.replace('_', ' ').upper() for name in gene_name_2]

    gene_co_occurrence_tuple_list = []
    for idx, co_occurence_times in enumerate(gene_co_occurrence):
        if co_occurence_times > 0:
            gene_co_occurrence_tuple_list.append((gene_name_1[idx], gene_name_2[idx], gene_co_occurrence[idx]))

    dict = {'gene_name_1': gene_name_1,
            'gene_name_2': gene_name_2,
            'co_occurrence': gene_co_occurrence}
    df_output = pd.DataFrame(dict)
    df_output.to_csv(output_folder + pmcid + '_co_occurrence_same_sent.csv', index=False)

    return gene_co_occurrence_tuple_list


def gene_co_occurrence_in_neighboring_sentence(sentence_list, gene_list, pmcid, output_folder):
    '''
    :param sentence_list: list, a list of sentences with each sentence as a string linked by '_'
    :param gene_list: list, a list of gene names
    :param csv_folder_path: str, a folder path to store gene co-occurrence csv file
    :param pmcid: str, pmcid, refer to the output gene co-occurrence csv file name
    :return: a dict, which contains gene co-occurrence information
                     gene_name_1, gene_name_2, co_occurrence, each as list
    '''

    occurrences = OrderedDict((name, OrderedDict((name, 0) for name in gene_list)) for name in gene_list)
    # A B C D E
    # A -> A B C D E
    # B -> A B C D E
    # Find the co-occurrences:
    neighbor = 0
    if neighbor >= len(sentence_list) - neighbor:
        tmp_gene_co_occur = {}
        three_sent = '_'.join(sentence_list)
        # print(three_sent)
        for gene in gene_list:
            # gene_replace = gene.replace(')', '\)').replace('(', '\(').replace(']', '\]').replace('[', '\[')
            single_gene_occur = [m.start() for m in re.finditer(gene, three_sent)]
            # ori_gene = gene_replace.replace('\\', '')
            tmp_gene_co_occur[gene] = len(single_gene_occur)
        # print(tmp_gene_co_occur)
        for gene_name, occur in tmp_gene_co_occur.items():
            for name in gene_list:
                occurrences[gene_name][name] += min(tmp_gene_co_occur[gene_name], tmp_gene_co_occur[name])
                occurrences[name][gene_name] += min(tmp_gene_co_occur[gene_name], tmp_gene_co_occur[name])
    else:
        for i in range(neighbor, len(sentence_list)-neighbor):
            tmp_gene_co_occur = {}
            three_sent = '_'.join(sentence_list[i-neighbor : i+neighbor+1])
            # print(three_sent)
            for gene in gene_list:
                # gene_replace = gene.replace(')', '\)').replace('(', '\(').replace(']', '\]').replace('[', '\[')
                single_gene_occur = [m.start() for m in re.finditer(gene, three_sent)]
                # ori_gene = gene_replace.replace('\\', '')
                tmp_gene_co_occur[gene] = len(single_gene_occur)
            # print(tmp_gene_co_occur)
            for gene_name, occur in tmp_gene_co_occur.items():
                for name in gene_list:
                    occurrences[gene_name][name] += min(tmp_gene_co_occur[gene_name], tmp_gene_co_occur[name])
                    occurrences[name][gene_name] += min(tmp_gene_co_occur[gene_name], tmp_gene_co_occur[name])

    # print(' ', ' '.join(occurrences.keys()))

    gene_name_1 = []
    gene_name_2 = []
    gene_co_occurrence = []
    for name_1, values in occurrences.items():
        for name_2, co_occurrence in values.items():
            if (co_occurrence != 0) and not (name_2 in gene_name_1 and name_1 in gene_name_2):
                gene_name_1.append(name_1)
                gene_name_2.append(name_2)
                gene_co_occurrence.append(co_occurrence)

    print(gene_name_1)
    print(gene_name_2)
    print(gene_co_occurrence)

    gene_name_1 = [name.replace('_', ' ').upper() for name in gene_name_1]
    gene_name_2 = [name.replace('_', ' ').upper() for name in gene_name_2]

    gene_co_occurrence_tuple_list = []
    gene_name_1_ = []
    gene_name_2_ = []
    gene_co_occurrence_ = []
    for idx, co_occurence_times in enumerate(gene_co_occurrence):
        if co_occurence_times > 0:
            gene_co_occurrence_tuple_list.append((gene_name_1[idx], gene_name_2[idx], gene_co_occurrence[idx]))
            gene_name_1_.append(gene_name_1[idx])
            gene_name_2_.append(gene_name_2[idx])
            gene_co_occurrence_.append(str(gene_co_occurrence[idx]))

    print(gene_name_1_)
    print(gene_name_2_)
    print(gene_co_occurrence_)


    dict = {'gene_name_1': gene_name_1_,
            'gene_name_2': gene_name_2_,
            'co_occurrence': gene_co_occurrence_}
    df_output = pd.DataFrame(dict)
    df_output.to_csv(output_folder + pmcid + '_co_occurrence_neighboring.csv', index=False)
    return gene_co_occurrence_tuple_list


def openie_match(full_text, gene_list, pmcid, output_folder):
    '''
    :param full_text: str, full text
    :param gene_set: set, contains all gene names from pathway figure
    :return: relation_tuples: a tuple list, confidence_score, gene_name_1, gene_name_2, relation
    '''

    # a tuple list -- confidence_score, gene_name_1, gene_name_2, relation
    relation_tuples = []
    sent_list = nltk.sent_tokenize(full_text)

    bi_gene = itertools.product(gene_list, gene_list)
    gene_pairs = [pair for pair in bi_gene if pair[0] != pair[1]]

    relations = []
    for i in range(10, len(sent_list)-10):
        ten_sent = '_'.join(sent_list[i - 10: i + 10 + 1])
        sent_file_name = './tmp.txt'
        with open(sent_file_name, "w") as text_file:
            text_file.write(ten_sent)

        process = Popen(['java', '-mx1g', '-cp', '*', 'edu.stanford.nlp.naturalli.OpenIE', sent_file_name], stdout=PIPE,
                        stderr=PIPE)
        stdout, stderr = process.communicate()
        relations += stdout.decode("utf-8").split('\n')

    confidence_score = []
    subject = []
    relation_entity = []
    object = []

    for relation in relations:
        for gene_pair in gene_pairs:
            if gene_pair[0] in relation and gene_pair[1] in relation:
                relation_list = list(relation.split('\t'))
                confidence_score.append(relation_list[0])
                subject.append(relation_list[1])
                relation_entity.append(relation_list[2])
                object.append(relation_list[3])
                # relation_tuples.append((relation_list[0], relation_list[1], relation_list[3], relation_list[2]))

    dict = {'confidence_score': confidence_score,
            'subject': subject,
            'relation': relation_entity,
            'object': object}
    df_output = pd.DataFrame(dict)
    df_output.to_csv(output_folder + pmcid + '_openie.csv', index=False)
    # return relation_tuples

#
# def dependency_parsing_match():
#

if __name__ == '__main__':
    # original data
    filename_pmcid_pair_file = 'data/figurename.xlsx'
    ground_truth_gene_names_file = 'data/finalized_genes.csv'
    ground_truth_gene_relations_file = 'data/finalized_relations.csv'
    # ocr_results_folder = 'data/preliminary_annotation_json/'
    whole_article_folder = 'data/full_text/'

    # output relations folder
    output_ocr_relation_folder = 'data/ocr_relations/'
    output_co_occur_same_folder = 'data/co_occur_same_sent/'
    output_co_occur_neigh_folder = 'data/full_text_one_text/'
    output_openie_folder = 'data/openie/'

    # figure name & corresponding pmcid xlsx
    df_filename_pmcid_pair = pd.read_excel(filename_pmcid_pair_file)
    cols = ['filename', 'PMCID']
    df_filename_pmcid_pair[cols] = df_filename_pmcid_pair[cols].apply(lambda x: x.str.strip())
    figure_names = list(df_filename_pmcid_pair['filename'])

    # ground truth csv -- gene names & relations
    df_gene_names_ground_truth = pd.read_csv(ground_truth_gene_names_file)
    df_gene_relations_ground_truth = pd.read_csv(ground_truth_gene_relations_file)

    # ground_truth_in_total = []

    for jpg_file in figure_names:
        # print(jpg_file)
        jpg_file_name = os.path.basename(jpg_file[:-4])
        # ocr_relations_file_name = ocr_results_folder + jpg_file_name + '_relation.json'
        pmcid = jpg_file_name.split('_')[0]

        # ground truth gene names
        df_ground_truth_gene_names = df_gene_names_ground_truth.loc[df_gene_names_ground_truth['fig_name'] == jpg_file]
        ground_truth_gene_names = list(set(df_ground_truth_gene_names['annotated_gene_name'].str.upper()))
        # ground_truth_in_total.append(str(len(ground_truth_gene_names)))

        # ground truth gene relations
        df_ground_truth_gene_relations = df_gene_relations_ground_truth.loc[df_gene_relations_ground_truth['fig_name'] == jpg_file]

        # # ocr results relations
        # startor_receptor_relation_tuple_list = parse_json(ocr_relations_file_name, pmcid, output_ocr_relation_folder)
        # # ocr_in_total.append(str(len(ocr_results_gene_names)))

        # whole article
        # whole_article = ''
        # with open(whole_article_folder + pmcid + '.txt', "r") as f:
        #     whole_article += f.read()
        #
        # # # gene co-occurence relations -- same & neighbor
        # processed_sent_list, processed_gene_list = preprocess_sent_list_and_gene_list(whole_article, ground_truth_gene_names)
        # # gene_co_occurrence_tuple_list = co_occurrence_same_sentence_match(processed_sent_list, ground_truth_gene_names, output_co_occur_same_folder)
        # gene_co_occurrence_tuple_list_3_sentences = gene_co_occurrence_in_neighboring_sentence(processed_sent_list, ground_truth_gene_names, pmcid, output_co_occur_neigh_folder)

        # openie extracted relations
        # openie_match(whole_article, ground_truth_gene_names, pmcid, output_openie_folder)
        # print(relation_tuples)
        print(",{}".format(len(df_ground_truth_gene_names)))
