import os
import json

if __name__ == "__main__":
    acutal_folder = r"E:\PythonDev\mathmodel\nlp_end_work\nlp_code\tie_data\hugo\parag_fig/"
    will_update_folder = r"E:\PythonDev\mathmodel\xml_help_full_text\tag_json\parag_fig/"
    for root, dirs, files in os.walk(acutal_folder):
        for file in files:
            with open(root + file, "r", encoding='utf-8') as r1:
                json_file_acutal = json.load(r1)
                only_key = ""
                for key in json_file_acutal.keys():
                    only_key = key
                    break
                with open(will_update_folder + file, "r", encoding='utf-8') as r2:
                    json_file_update = json.load(r2)
                    flag = False
                    for (key, value) in json_file_update.items():
                        if key == only_key:
                            flag = True
                            with open("process_para/" + file.split(".")[0] + ".txt", "w", encoding='utf-8') as w :
                                w.write(value)
                            break
                    if not flag:
                        print("Debug: {}".format(file))