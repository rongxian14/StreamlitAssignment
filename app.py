import streamlit as st
import pandas as pd
import networkx as nx
import requests
from Bio import Entrez
from Bio import SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import matplotlib.pyplot as plt

#1. Functions to get protein characteristics
# To retrieve protein sequence given Uniprot ID
def retSeq(id):
    handle = Entrez.efetch(db='protein', id=id, rettype = 'fasta', retmode = 'text')
    output = SeqIO.read(handle, 'fasta')
    return output

# To get the protein length
def protLength(seq):
    return len(str(seq))

# To get the protein weight
def protWeight(seq):
    protein_analysis = ProteinAnalysis(str(seq))
    protWeight = round(protein_analysis.molecular_weight(),2)
    return protWeight

def get_protein_name(protName):
    # Split the string by "|"
    parts = protName.split("|")
    # Extract the word between "|" and "_"
    if len(parts) >= 3:
        protName = parts[2].split("_")[0]
        return protName

def get_protein_characteristics(protein_data):
    # To get the protein name
    protName = get_protein_name(protein_data.name)
   
    # To get the protein length
    protLength = len(str(protein_data.seq))

    # The Protein Analysis module from Bio.SeqUtils.ProtParam is used...
    protein_analysis = ProteinAnalysis(str(protein_data.seq))
    # To get the protein weight
    protWeight = round(protein_analysis.molecular_weight(),2)
    # To get the protein isoelectric point
    protIsoelectricPoint = round(protein_analysis.isoelectric_point(),2)
    # To get the amino acid composition
    aa_composition = protein_analysis.get_amino_acids_percent()
    return protName, protLength, protWeight,protIsoelectricPoint, aa_composition
   
#2. Functions to plot PPI graph
def retrieve_ppi(protein_name):
  string_url = "https://string-db.org/api/json/network"
  params = {
      "identifiers":protein_name,
      "species":9606
  }
  response = requests.get(string_url, params=params)
  data = response.json() #parse to json
  network_df = pd.json_normalize(data)
  return network_df

def get_characteristics(ppi_graph):
    st.write("Number of nodes:", ppi_graph.number_of_nodes())
    st.write("Number of edges:", ppi_graph.number_of_nodes())
    st.write("Number of interactions of each nodes:", ppi_graph.degree())
    degree_central = nx.degree_centrality(ppi_graph)
    top_5_proteins = sorted(degree_central.items(), key=lambda x:-x[1])[:5]
    high_centrality_proteins = [node for node, centrality in top_5_proteins]
    st.write(f"Top 5 proteins: {high_centrality_proteins}")
    return high_centrality_proteins

def visualize_ppi(protein_name):
    df = retrieve_ppi(protein_name)
    ppi_graph = nx.from_pandas_edgelist(df, "preferredName_A", "preferredName_B")
    graph = nx.spring_layout(ppi_graph)
    high_centrality_proteins = get_characteristics(ppi_graph)
    nx.draw(ppi_graph, graph, with_labels=True, node_size=500, node_color='pink', edge_color='blue', font_size=8)
    nx.draw_networkx_nodes(ppi_graph, graph, nodelist=high_centrality_proteins, node_color='orange') #the top 5 will be orange
    plt.title("Protein-Protein Interaction Network")
    plt.axis("off")  # Hide axis
    st.pyplot()  # Display the graph in Streamlit
    st.set_option('deprecation.showPyplotGlobalUse', False)

#3. Functions for analysis of protein sequences
def folded_or_unfolded(protein_seq):

    # Kyte-Doolittle hydrophobicity scale
    aa_hydrophobicity = {'I': 4.50, 'V': 4.20, 'L': 3.80, 'F': 2.80, 'C': 2.50, 'M': 1.90, 'A': 1.80, 'G': -0.40, 'T': -0.70, 'S': -0.80,
                         'W': -0.90,'Y': -1.30, 'P': -1.60, 'H': -3.20, 'E': -3.50, 'N': -3.50, 'Q': -3.50, 'D': -3.50, 'K': -3.90, 'R': -4.50}

    # Version 2 (Brute-Force) of calculating the mean hydrophobicity values:
    mean_hydrophobicity = 2*(aa_hydrophobicity[protein_seq[0]] + aa_hydrophobicity[protein_seq[1]] + aa_hydrophobicity[protein_seq[-1]] + aa_hydrophobicity[protein_seq[-2]] + ((len(protein_seq)-4)/5))/(len(protein_seq))

    # Compares the mean hydrophobicity value calculated with the experimentally-derived hydrophobicity cutoff values of natively folded and unfolded proteins
    if (mean_hydrophobicity >= 0.34 and mean_hydrophobicity <= 0.44):
        return (f"Since the mean hydrophobicity of the protein is {mean_hydrophobicity} which is >=0.34 and <=0.44. The protein is natively unfolded.")
    elif (mean_hydrophobicity >= 0.45 and mean_hydrophobicity <= 0.51):
        return (f"Since the mean hydrophobicity of the protein is {mean_hydrophobicity} which is >=0.45 and <=0.51. The protein is natively folded.")
    else:
        return ('Inconclusive')

def plot_hydrophobicity_profile(protein_seq, smoothing_window=1):

    # This function serves to plot a hydrophobicity profile for a protein sequence
    # X axis has the amino acids, sequentially arranged
    # Y axis has the hydrophobic value for that particular amino acid
    # These dots are then connected to form the hydrophobic profile

    import matplotlib.pyplot as plt
    aa_hydrophobicity = {'I': 4.50, 'V': 4.20, 'L': 3.80, 'F': 2.80, 'C': 2.50, 'M': 1.90, 'A': 1.80, 'G': -0.40, 'T': -0.70, 'S': -0.80,
                         'W': -0.90,'Y': -1.30, 'P': -1.60, 'H': -3.20, 'E': -3.50, 'N': -3.50, 'Q': -3.50, 'D': -3.50, 'K': -3.90, 'R': -4.50}
    plt.title('Hydrophobicity Profile')
    plt.xlabel('Amino acid sequence')

    x_values = []
    y_values = []
    if smoothing_window == 1:
        plt.ylabel('Hydrophobicity index')
        for aa in range(len(protein_seq)):
            x_values.append(aa)
            y_values.append(aa_hydrophobicity[protein_seq[aa]])
        plt.plot(x_values, y_values)

    # Smoothing window plots the average hydrophobicity value of the amino acids within the odd window

    elif smoothing_window > 1:
        for aa in range(len(protein_seq)):
            y_values.append(aa_hydrophobicity[protein_seq[aa]])
        half_window = smoothing_window//2
        plt.ylabel('Hydrophobicity index (window = %d)' % smoothing_window)
        y_average = []
        for i in range(half_window, len(protein_seq)-half_window):
            average = 0.0
            x_values.append(i)
            for j in range(-half_window, half_window+1):
                average += y_values[i+j]
            y_average.append(average/smoothing_window)
        plt.plot(x_values, y_average)
    st.pyplot(plt)
    

# Streamlit App
st.title('Protein Analysis App')

# Sidebar
input_type = st.sidebar.radio('Input Type', ['Uniprot ID', 'Protein Sequence'])

if input_type == 'Uniprot ID':
    uniprot_id = st.sidebar.text_input('Enter Uniprot ID')
    if uniprot_id == '':
        st.error('Please enter a Uniprot ID, such as P00533')
    
    if st.sidebar.button('Fetch Protein Data'):
        # Retrieve protein data using Uniprot ID
        protein_data = retSeq(uniprot_id)
                
        # Display protein characteristic
        st.subheader('Protein Characteristics')
        df = pd.DataFrame(
            {
                "Characteristics": ["Name", "Length", "Weight(Daltons)", "Isoelectric Point(pI)","Amino Acid Composition"],
                "Description": [get_protein_characteristics(protein_data)[0], 
                                get_protein_characteristics(protein_data)[1], 
                                get_protein_characteristics(protein_data)[2], 
                                get_protein_characteristics(protein_data)[3],
                                get_protein_characteristics(protein_data)[4],
                ],
            }
        )
        st.dataframe(df, hide_index=True)
        st.divider()
        
        # Display protein-protein interaction network
        st.subheader('Protein-Protein Interaction Network')
        
        protein_name = str(get_protein_name(protein_data.name))
        st.write(f"Visualizing PPI network for protein {protein_name}...")
        visualize_ppi(protein_name)
        st.set_option('deprecation.showPyplotGlobalUse', False)

elif input_type == 'Protein Sequence':
    protein_sequence = st.sidebar.text_area('Enter Protein Sequence')
    if protein_sequence == '':
        st.error('Please enter a protein sequence, such as MGLALLLLAVLGATATEQSSLCAQVHGQICCTHGGCSGTSIQCLPCHSQCAAGCTGPKIPSIATGMVGALLLLLVVALGIGLFMRRRHIVRKRTLRRLLQERELVEPLTPSGEAPNQALLRILKETEAYVMASVDNPHVCRLLGICLTSTVQLITQLMPFGCLLDYVREHKDNIGSQYCVLSHNVVKHERCEQTKNGQGCKYVNENKSWHCPPICYVALLNKDGKVCINADGRYEKMSKCGAPDCVKTQVCANRCSGRCWGGCVRQCIQFAYELAECCQPELLPLGQKANKEALQKCEE')
    if st.sidebar.button('Analyze Sequence'):
    # To anayse sequences:
        # 1. Categorization (natively folded vs. natively unfolded) according to mean hydrophobicity
        st.subheader(f"Calculating the mean hydrophobicity of the protein...")
        st.write(folded_or_unfolded(protein_sequence))
        # 2. Secondary structure prediction
        st.subheader(f"Plotting the hydrophobicity profile of the protein...")
        plot_hydrophobicity_profile(protein_sequence.upper(), 5)
    
