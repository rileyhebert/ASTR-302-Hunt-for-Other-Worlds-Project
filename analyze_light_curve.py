# -*- coding: utf-8 -*-
"""
Created on Wed Nov 18 10:32:10 2026

@author: Riley
"""

import csv
import matplotlib.pyplot as plt
import numpy as np


def get_data(output_file):
    list_of_dicts = []
    with open(output_file,"r",newline = "") as file:
        reader = csv.DictReader(file)
        for row in reader:
            list_of_dicts.append(row)
    return list_of_dicts

def graph_fluxes_with_errors(data_list):
    x = np.arange(0,len(data_list),1)
    y = list(data_list[0].values())
    y_err = list(data_list[1].values())
    print(data_list[2].values())
    plt.errorbar(x, y,yerr=y_err)
    plt.xlabel("Image Number")
    plt.ylabel("Flux (counts)")
    plt.title("Flux of the Star Over Time")
    plt.show()
    
def main():
    output_file = input("Name of output file:\n")
    data = get_data(output_file)
    graph_fluxes_with_errors(data)
    
main()