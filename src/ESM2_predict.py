#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 12 15:15:40 2023

@author: paceramateos
"""

import pandas as pd
pd.set_option('display.max_columns', None)
import torch
import esm
import pandas as pd
import requests
import numpy as np
import plotly.graph_objects as go
import plotly.offline as offline
import argparse

def get_uniprot_sequence(uniprot_id):
    url = f'https://www.uniprot.org/uniprot/{uniprot_id}.fasta'
    response = requests.get(url)

    if response.status_code != 200:
        print(f"Failed to retrieve data for {uniprot_id}. Status code: {response.status_code}")
        return None

    fasta_data = response.text
    sequence = "".join(fasta_data.split("\n")[1:])
    return sequence


def main():    
    parser = argparse.ArgumentParser(description="Script to process UniProt ID and generate outputs")
    parser.add_argument("-u", "--uniprot_id", required=True, help="UniProt ID to be processed")
    parser.add_argument("-f", "--model_folder", required=True, help="Folder path where the mode")
    parser.add_argument("--plot_llr", default=True, help="Flag to indicate if the log likelihood ratio should be plotted")
    parser.add_argument("--plot_heatmap", default=True, help="Flag to indicate if the heatmap should be plotted")
    
    args = parser.parse_args()

    uniprot_id = args.uniprot_id
    model_folder = args.model_folder
    plot_llr = args.plot_llr
    plot_heatmap = args.plot_heatmap
    
    sequence_aa = get_uniprot_sequence(uniprot_id)
    
    model_path = model_folder
    
    # Load ESM-2
    model, alphabet = esm.pretrained.load_model_and_alphabet_local(model_path +'/esm2_t33_650M_UR50D.pt')
    
    batch_converter = alphabet.get_batch_converter()
    
    # Load ESM-2 model
    model.eval()  # disables dropout for deterministic results
    
    data = list(zip([uniprot_id], [sequence_aa]))
    
    # Convert input sequences to tokens
    batch_labels, batch_strs, batch_tokens = batch_converter(data)
    batch_lens = (batch_tokens != alphabet.padding_idx).sum(1)
    
    # Stored results
    results = {}
    for i, batch in enumerate(batch_tokens):                 
        results[batch_labels[i]] = model(batch.unsqueeze(0), repr_layers=[33], return_contacts=True)
    
    
    # converts logits in probabilities
    for key in enumerate(results.keys()):
    
        probabilities = torch.nn.functional.softmax(results[key[1]]['logits'], dim=-1)
        
        log_like_ratio = []
        for pos in range(1,int(batch_lens[key[0]])-1):
            position_in_sequence = []
            for aa in enumerate(alphabet.standard_toks[:-7]):
                position_in_sequence.append(float(np.log(probabilities[0][pos][alphabet.get_idx(aa[1])].squeeze(0).detach()/probabilities[0][pos][alphabet.get_idx(data[key[0]][1][pos-1].upper())].squeeze(0).detach())))
            log_like_ratio += [position_in_sequence]
        
        # Plot heatmap
        
        # add a flag if want to plot the log likelihood ratio
        if plot_llr:
            fig = go.Figure(go.Heatmap(z=np.array(log_like_ratio).T, colorscale='RdBu'))
            
            yticks_values = list(range(len(alphabet.standard_toks[:-7])))  # Generate ytick values based on the number of labels
            yticks_labels = ['L', 'A', 'G', 'V', 'S', 'E', 'R', 'T', 'I', 'D', 'P', 'K', 'Q', 'N', 'F', 'Y', 'M', 'H', 'W', 'C']
            
            # Update the y-axis settings with custom yticks
            fig.update_yaxes(tickvals=yticks_values, ticktext=yticks_labels)
            
            fig.update_layout(title=key[1])
            
            # Show the heatmap
            offline.plot(fig)
        
    
    # Plot the heatmap
    for key in enumerate(results.keys()):
            
        # Assuming you have already computed the heatmap data
        heatmap_data = results[key[1]]['contacts'].squeeze(0).detach()[: len(batch_strs[key[0]]), :len(batch_strs[key[0]])]
        if plot_heatmap:
    
            # Create the heatmap using Plotly
            fig = go.Figure(data=go.Heatmap(z=heatmap_data, colorscale='Reds'))
            
            # Set the title
            fig.update_layout(title=key[1])
            
            # Display the plot in the browser
            offline.plot(fig)
    
    
    # Save the data to files
    with open(f"{uniprot_id}_heatmap_data.txt", "w") as heatmap_file:
       heatmap_file.write(str(heatmap_data))

        
    with open(f"{uniprot_id}_log_likelihood_ratio.txt", "w") as llr_file:
       llr_file.write(str(log_like_ratio))
        
if __name__ == "__main__":
    main()