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
from torch import nn
import pandas as pd
import requests
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import plotly.express as px
import re
from scipy.stats import gmean
from scipy.stats import pearsonr
import ast
import plotly.graph_objects as go
import plotly.offline as offline
import pickle


def get_uniprot_sequence(uniprot_id):
    url = f'https://www.uniprot.org/uniprot/{uniprot_id}.fasta'
    response = requests.get(url)

    if response.status_code != 200:
        print(f"Failed to retrieve data for {uniprot_id}. Status code: {response.status_code}")
        return None

    fasta_data = response.text
    sequence = "".join(fasta_data.split("\n")[1:])
    return sequence


sequence_aa = get_uniprot_sequence('Q9UHP9')


# Load ESM-2
model, alphabet = esm.pretrained.load_model_and_alphabet_local('/g/data/th81/mutation_prediction/models/esm2_t33_650M_UR50D.pt')

batch_converter = alphabet.get_batch_converter()

# Load ESM-2 model
model.eval()  # disables dropout for deterministic results

data = list(zip(['SMARCA2', 'SMARCA2_criptic'], [SMARCA2, SMARCA2_criptic]))

# Convert input sequences to tokens
batch_labels, batch_strs, batch_tokens = batch_converter(data)
batch_lens = (batch_tokens != alphabet.padding_idx).sum(1)









