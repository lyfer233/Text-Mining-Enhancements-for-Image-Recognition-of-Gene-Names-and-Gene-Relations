# -*- coding: utf-8 -*-
# @Time    : 2021/12/13 12:10

import os
import xml.etree.ElementTree as ET
from xml.etree import ElementTree
from Bio import Entrez
import os
import json
import pandas as pd
import re
from collections import defaultdict

def pcmid_to_xml_file(PMCIDs, save_xml_file):
    for i, pmc_id in enumerate(PMCIDs):
        file_path = "{}/{}.xml".format(save_xml_file, pmc_id)
        if not os.path.isfile(file_path):
            net_handle = Entrez.efetch(db="pmc", id=pmc_id, rettype="medline", retmode="xml")
            # net_handle = Entrez.efetch(db="pmc", id="PMC6681686", rettype="medline", retmode="text")
            out_handle = open(file_path, "wb")
            out_handle.write(net_handle.read())
            out_handle.close()
            net_handle.close()
            print("{}. {} Saved".format(i, pmc_id))

def parse_paragraph(dir_path, file_name):
    file = os.path.join(dir_path, file_name)
    tree = ET.parse(file)
    root = tree.getroot()
    infor = {}
    for each_tag in root.iter('sec'):
        # for each_tag in object_tag:
        label = each_tag.find("title").text
        infor[label] = {}
        for i, para in enumerate(each_tag.findall("p")):
            temp = []
            for item in para.itertext():
                if item != '\n':
                    temp.append(item)
            infor[label]['p'+str(i)] = " ".join(temp)
    # output = {}
    # for key, values in infor.items():
    #     output[key] = " ".join(values)
    return infor

def parse_full(dir_path, file_name):
    file = os.path.join(dir_path, file_name)
    tree = ET.parse(file)
    root = tree.getroot()
    infor = ''
    for para in root.iter('p'):
        for item in para.itertext():
            if item.strip():
                infor += item.strip()
    for each_tag in root.iter('sec'):
        # for each_tag in object_tag:
        try:
            label = each_tag.find("title").text
            infor += label.strip()
        except:
            continue
        for para in each_tag.findall("p"):
            for item in para.itertext():
                if item.strip():
                    infor += item.strip()
    return infor


def parse_caption(dir_path, file_name):
    file = os.path.join(dir_path, file_name)
    tree = ET.parse(file)
    root = tree.getroot()
    infor = {}
    for each_tag in root.iter('fig'):
    # for each_tag in object_tag:
        label = each_tag.find("label").text
        infor[label] = []
        for caption in each_tag.findall("caption"):
            for para in caption.findall("p"):
                for item in para.itertext():
                    if item.strip():
                        infor[label].append(item.strip())
    output = {}
    for key, values in infor.items():
        output[key] = " ".join(values)
    return output

def parse_fig_paragraph(dir_path, file_name):
    file = os.path.join(dir_path, file_name)
    tree = ET.parse(file)
    root = tree.getroot()
    infor = defaultdict(list)
    pattern = re.compile(r'\d+')

    for para in root.iter('p'):
        temp = []
        for item in para.itertext():
            if item.strip():
                temp.append(item.strip())
        if 'fig' in " ".join(temp) or "Fig" in " ".join(temp):
            m = re.findall(r"fig.+\d", " ".join(temp), re.I)
            for _ in m:
                # items = pattern.findall(_.split())
                items = _.split()
                fig_num = []
                for i, i_ in enumerate(items):
                    if 'fig' in i_ or 'Fig' in i_:
                        fig_num.append(items[i + 1])
                    try:
                        if ('fig' in i_ or 'Fig' in i_) and items[i + 2] == 'and':
                            fig_num.append(items[i + 3])
                    except:
                        continue
                for fn in fig_num:
                    infor['fig ' + fn].append(" ".join(temp))
                    # infor['fig ' + fn].append(temp)

    for each_tag in root.iter('sec'):
        # for each_tag in object_tag:
        for i, para in enumerate(each_tag.findall("p")):
            temp = []
            for item in para.itertext():
                if item.strip():
                    temp.append(item.strip())
            if 'fig' in " ".join(temp) or "Fig" in " ".join(temp):
                m = re.findall(r"fig.+\d+", " ".join(temp), re.I)
                for _ in m:
                    # items = pattern.findall(_.split())
                    items = _.split()
                    fig_num = []
                    for i, i_ in enumerate(items):
                        if 'fig' in i_ or 'Fig' in i_:
                            fig_num.append(items[i+1])
                        try:
                            if ('fig' in i_ or 'Fig' in i_) and items[i+2] == 'and':
                                fig_num.append(items[i + 3])
                        except:
                            continue
                    for fn in fig_num:
                        infor['fig '+fn].append(" ".join(temp))
    output = {}
    for key, values in infor.items():
        output[key] = " ".join(set(values))
    print(output)
    return output

if __name__ == '__main__':

    # file = '45_all_text.csv'
    # csvPD = pd.read_csv(file)
    # PMCIDs = {name.split('__')[0] for name in csvPD['fig_name']}
    # # PMCIDs = ['PMC2172299', 'PMC2659383', 'PMC3068308']
    # # PMCIDs = ['PMC3789177', 'PMC4061041']
    # print("PMCID num: {}".format(len(PMCIDs)))
    # Entrez.email = "qujing579@163.com"
    # save_xml_file = './xml_file'
    # if not os.path.isdir(save_xml_file):
    #     os.mkdir(save_xml_file)
    # pcmid_to_xml_file(PMCIDs, save_xml_file)
    #
    save_xml_file = './xml_file'
    file_names = [dir for dir in os.listdir(save_xml_file)]

    save_parsed_floder = "./paragraph"
    if not os.path.isdir(save_parsed_floder):
        os.mkdir(save_parsed_floder)
    n = 0
    tag = 'parag_fig'
    for name in file_names:
        # print(name)
        tag_dict = {}
        if tag == 'caption':
            tag_dict = parse_caption(save_xml_file, name)
        elif tag == 'parag_fig':
            tag_dict = parse_fig_paragraph(save_xml_file, name)
        elif tag =='full':
            tag_dict = parse_full(save_xml_file, name)

        save_tag_floder = './tag_json/' + tag
        if not os.path.isdir(save_tag_floder):
            os.mkdir(save_tag_floder)
        if len(tag_dict):
            n += 1
            with open("{}/{}.json".format(save_tag_floder, name.split('.')[0]), 'w', encoding='utf-8') as jf:
                json.dump(tag_dict, jf, ensure_ascii=False, indent=4)
    print(n)