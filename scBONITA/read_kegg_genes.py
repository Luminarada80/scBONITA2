import requests
import logging

def get_gene_names_from_kegg_pathway(pathway_id):
    base_url = "http://rest.kegg.jp"
    get_pathway_url = f"{base_url}/get/{pathway_id}"
    response = requests.get(get_pathway_url)
    
    if response.status_code != 200:
        print("Error fetching pathway data")
        return []
    
    pathway_data = response.text
    gene_names = []

    start = False
    for line in pathway_data.split('\n'):
        if line.startswith("GENE"):
            start = True
        
        if start == True:
            if ";" in line:
                csv_format = line.replace("GENE","\t").replace("        ", "\t").replace("    ", "\t").replace("  ", "\t").replace("; ", "\t")
                gene_info = csv_format.split('\t')[3]
                gene_names.append(gene_info)

            else:
                start = False
        
    return gene_names

if __name__ == '__main__':

    # Set the logging level for output
    logging.basicConfig(format='%(message)s', level=logging.INFO)

    kegg_pathway = input("Which KEGG pathway are you interested in?: ")
    gene_names = get_gene_names_from_kegg_pathway(kegg_pathway)
    print(f'\nGenes:')
    print(gene_names)

    print(f'\nWriting to file "{kegg_pathway}.csv"...')

    with open(f'{kegg_pathway}.csv', 'w') as outfile:
        outfile.write(f'Pathway: {kegg_pathway}\n')
        for gene in gene_names:
            outfile.write(f'{gene}\n')