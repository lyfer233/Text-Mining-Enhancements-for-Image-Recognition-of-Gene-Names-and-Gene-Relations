import os
import json
if __name__ == "__main__":
    json_file_folder = "./tag_json/"
    full_text_folder = "../yijie/data/full_text/"
    save_folder = "../yijie/data/figure/"
    for root, dirs, files in os.walk(json_file_folder):
        for file in files:
            PMCID = str(file).split('_')[0]
            with open(json_file_folder + file, "r", encoding='utf-8') as r1:
                json_file = json.load(r1)
                figure_file = ""
                for tmp in json_file.values():
                    figure_file = tmp
                    break
                with open(full_text_folder + PMCID + ".txt", "r", encoding='utf-8') as r2:
                    full_text_file = r2.read()
                    # print(PMCID, end="")
                    if figure_file in full_text_file:
                        with open(save_folder + PMCID + ".txt", "w", encoding='utf-8') as w:
                            print(PMCID)
                            w.write(figure_file)
                    # else:
                    #     print(" no", end="  ")
                    #     print(figure_file)