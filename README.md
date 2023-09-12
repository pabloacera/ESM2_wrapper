# ESM2_wrapper

### ESM2 installation
pip install fair-esm  # latest release

mkdir models
cd models

**1. download the ESM2 model**
```bash
wget https://dl.fbaipublicfiles.com/fair-esm/models/esm2_t33_650M_UR50D.pt
```

**2. download the regression model for the contact map**
```bash
wget https://dl.fbaipublicfiles.com/fair-esm/regression/esm2_t33_650M_UR50D-contact-regression.pt
```

**3. Install Dependencies**

```bash
pip install -r requirements.txt
```

**4. Information about how to use the script**

```bash
python ESM2_predict.py -h
usage: ESM2_predict.py [-h] -u UNIPROT_ID -f MODEL_FOLDER [--plot_llr PLOT_LLR] [--plot_heatmap PLOT_HEATMAP]

Script to process UniProt ID and generate outputs

optional arguments:
  -h, --help            show this help message and exit
  -u UNIPROT_ID, --uniprot_id UNIPROT_ID
                        UniProt ID to be processed
  -f MODEL_FOLDER, --model_folder MODEL_FOLDER
                        Folder path where the mode
  --plot_llr PLOT_LLR   Flag to indicate if the log likelihood ratio should be plotted
  --plot_heatmap PLOT_HEATMAP
                        Flag to indicate if the heatmap should be plotted
```

## Outputs

This scripts uses as input a UNIPROT_ID and will predict the likelihood ratio for each position for each aminoacid and also the residues interaction heatmap. 
It will generate two files called {UNIPROT_ID}_heatmap_data.txt and {UNIPROT_ID}_log_likelihood_ratio.txt containing the data from the heatmap and log likelihood ration respectivelly. 
By default will generate plotly plots for each of the two items.



