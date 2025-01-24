
"""
main_analysis.py
================

Main script for analyzing avian influenza mutation patterns using NCBI data.
"""

import os
from Bio import Entrez, SeqIO
import pandas as pd
from sklearn.cluster import KMeans
import matplotlib.pyplot as plt
import json

credentials = json.load("../credentials.json")
# Set up email for NCBI access
Entrez.email = credentials['email']

def fetch_ncbi_metadata(query, max_results=100):
    """
    Fetch metadata from NCBI using Biopython.

    Parameters:
        query (str): Search query for NCBI.
        max_results (int): Maximum number of results to fetch.

    Returns:
        str: GenBank records as a string.
    """
    handle = Entrez.esearch(db="nucleotide", term=query, retmax=max_results)
    record = Entrez.read(handle)
    handle.close()
    
    ids = record["IdList"]
    handle = Entrez.efetch(db="nucleotide", id=ids, rettype="gb", retmode="text")
    data = handle.read()
    handle.close()
    return data

def parse_metadata(genbank_data):
    """
    Parse metadata from GenBank records into a DataFrame.

    Parameters:
        genbank_data (str): GenBank records as a string.

    Returns:
        pd.DataFrame: Metadata as a DataFrame.
    """
    records = SeqIO.parse(genbank_data.splitlines(), "genbank")
    data = []
    for record in records:
        annotations = {
            "accession": record.id,
            "organism": record.annotations.get("organism", ""),
            "date": record.annotations.get("date", ""),
            "host": record.annotations.get("host", "unknown"),
        }
        data.append(annotations)
    return pd.DataFrame(data)

def cluster_metadata(df):
    """
    Perform KMeans clustering on metadata.

    Parameters:
        df (pd.DataFrame): Metadata DataFrame.

    Returns:
        pd.DataFrame: DataFrame with cluster assignments.
    """
    df["host_encoded"] = pd.factorize(df["host"])[0]
    kmeans = KMeans(n_clusters=3, random_state=42)
    df["cluster"] = kmeans.fit_predict(df[["host_encoded"]])
    return df

def visualize_clusters(df):
    """
    Visualize clusters using matplotlib.

    Parameters:
        df (pd.DataFrame): DataFrame with cluster assignments.
    """
    plt.scatter(df["host_encoded"], df.index, c=df["cluster"], cmap="viridis")
    plt.title("Clusters of Bird Flu Strains by Host")
    plt.xlabel("Host (Encoded)")
    plt.ylabel("Strain Index")
    plt.show()

if __name__ == "__main__":
    # Example query: H5N1 strains
    query = '"H5N1"[Organism] AND "avian"[Host]'
    print("Fetching metadata...")
    metadata = fetch_ncbi_metadata(query)
    
    print("Parsing metadata...")
    df = parse_metadata(metadata)
    
    print("Clustering metadata...")
    df = cluster_metadata(df)
    
    print("Visualizing clusters...")
    visualize_clusters(df)
