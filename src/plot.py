import os
import json
import base64
import subprocess
import pandas as pd

def plot_fusions(
        gene_list,#type:ignore
        pdf_output#type:ignore
        ):
    message = ""
    json_handle = open("gene_to_coordinates.json", 'r')
    gene_dict = json.load(json_handle)
    increase = lambda x : int(x) + 10**(len(str(x))-2)#type:ignore
    decrease = lambda x : int(x) - 10**(len(str(x))-2)#type:ignore
    left_data = {"chr":[],"start":[],"end":[],"gene":[]}#type:ignore
    right_data = {"chr":[],"start":[],"end":[],"gene":[]}#type:ignore
    for genes in gene_list:#type:ignore
        if len(genes)==2:#type:ignore
            left_gene = genes[0]#type:ignore
            right_gene = genes[1]#type:ignore
            if left_gene in gene_dict.keys() and right_gene in gene_dict.keys():
                l_chromosome = gene_dict[left_gene][1]
                l_start = gene_dict[left_gene][2]
                l_end = gene_dict[left_gene][3]
                r_chromosome = gene_dict[right_gene][1]
                r_start = gene_dict[right_gene][2]
                r_end = gene_dict[right_gene][3]
                if l_chromosome != "MT" and r_chromosome != "MT":
                    left_data["chr"].append(f"chr{l_chromosome}")#type:ignore
                    left_data["start"].append(decrease(l_start))#type:ignore
                    left_data["end"].append(increase(l_end))#type:ignore
                    left_data["gene"].append(left_gene)#type:ignore
                    right_data["chr"].append(f"chr{r_chromosome}")#type:ignore
                    right_data["start"].append(decrease(r_start))#type:ignore
                    right_data["end"].append(increase(r_end))#type:ignore
                    right_data["gene"].append(right_gene)#type:ignore
                else:
                    message+=f"Left gene: {left_gene} or right gene: {right_gene} present on mitochondrial chromosome!\n"
                    #print(f"Left gene: {left_gene} or right gene: {right_gene} present on mitochondrial chromosome!")
            else:
                message+=f"Left gene: {left_gene} or right gene: {right_gene} not present in gtf!\n"
                #print(f"Left gene: {left_gene} or right gene: {right_gene} not present in gtf !")
    left_coordinates_dataframe = pd.DataFrame(left_data)
    right_coordinates_dataframe = pd.DataFrame(right_data)
    if not left_coordinates_dataframe.empty and not right_coordinates_dataframe.empty:
        left_coordinates_dataframe.to_csv("left.bed", index=False,sep="\t")
        right_coordinates_dataframe.to_csv("right.bed", index=False, sep="\t")
        command = [#type:ignore
                "Rscript",
                "plot.R",
                "-l",
                "left.bed",
                "-r",
                "right.bed",
                "-o",
                pdf_output
                ]
        run_status = subprocess.run(command)#type:ignore
        if os.path.exists(pdf_output):#type:ignore
            with open(pdf_output, 'rb') as file:#type:ignore
                pdf_data = file.read()
            base64_data = base64.b64encode(pdf_data).decode('utf-8')
            return {
            "data": base64_data,
            "message": message,
            "error_code":None,
            "success":True
            }
        else:
            return {
            "data": None,
            "message": message,
            "error_code":run_status,
            "success":False
            }
    else:
        #print("No data in left and right coordinate files !")
        message+="No data in left and right coordinate files!\n"
    return {
        "data": None,
        "message": message,
        "error_code":"no input provided!",
        "success":False
    }

# genes_fusion = []
# fusions = "20900134532_fusions.tsv"
# for line in open(fusions, 'r').readlines():
#     if "#FusionName" not in line:
#         l = [l.strip() for l in line.split()]
#         k=l[0].split("--")
#         genes_fusion.append([k[0], k[1]])#type:ignore

# #print(genes_fusion)
# a=plot_fusions(genes_fusion, pdf_output = "2new_output.pdf")
# print(a["message"])


    
