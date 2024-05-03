import pandas as pd
import os

def import_metadata(metadata_filepath, filtered_metadata_filepath):
    cell_names = [""]
    with open(metadata_filepath, 'r') as metadata_file, open(filtered_metadata_filepath, 'w') as filtered_metadata_file:
        for line in metadata_file:
            line = line.strip().split(',')
            name, group, cluster_name = line[0], line[6], line[8]
            if 'T cells CD8 - 1' in cluster_name:
                filtered_metadata_file.write(f'{name},{group}\n')
                cell_names.append(name)

    return cell_names

def clean_and_process_file(filepath, cell_names):
    """ Load, clean, and return the cleaned dataframe, filtering by cell names. """
    df = pd.read_csv(filepath)
    # # Handle an unnamed first column, assume it's the gene names
    # if df.columns[0] == '':
    #     df.rename(columns={ df.columns[0]: "Gene" }, inplace=True)
    # elif df.columns[0].strip() != 'Gene':
    #     # If the first column is not named "Gene" and is not empty, rename it
    #     df.rename(columns={ df.columns[0]: "Gene" }, inplace=True)

    # cols_to_keep = ['Gene'] + [col for col in df.columns[1:] if col in cell_names]  # Exclude the first column from the cell name check

    df = df.loc[:, df.columns.isin(cell_names)]  # Keep only columns listed in cell_names

    print(df.head())
    return df

def read_and_index_file(filepath, cell_names):
    """Read a file and set the first column as the index."""
    df = pd.read_csv(filepath, index_col=0)
    df = df.loc[:, df.columns.isin(cell_names)]
    return df

directory = '../../kazer_data/'
metadata_file = '../../kazer_data/merged_metadata.csv'
filtered_metadata_file = '../../kazer_data/filtered_metadata.csv'
output_file = '../../kazer_data/merged_data.csv'

dataframes = []

cell_names = import_metadata(metadata_file, filtered_metadata_file)


# Read each file, set the first column as index (assuming it contains gene names)
for filename in os.listdir(directory):
    if filename.startswith("GSM") and filename.endswith("expression_data.csv"):
        filepath = os.path.join(directory, filename)
        try:
            df = read_and_index_file(filepath, cell_names)
            dataframes.append(df)
        except Exception as e:
            print(f"Failed to process {filename}: {e}")

# Merge all dataframes along columns (side by side), aligning on the index (gene names)
if dataframes:
    merged_df = pd.concat(dataframes, axis=1, join='outer')
    merged_df.fillna(0, inplace=True)

    # Counting non-zero entries for each column
    column_gene_count = merged_df.apply(lambda x: (x > 0).sum(), axis=0)
    # Counting non-zero entries for each row
    row_gene_count = merged_df.apply(lambda x: (x > 0).sum(), axis=1)

    # Dropping columns where non-zero count is less than 200
    cols_to_drop = [col for col, count in column_gene_count.items() if count < 200]
    merged_df.drop(cols_to_drop, axis=1, inplace=True)

    # Dropping rows where non-zero count is less than 3
    rows_to_drop = [row for row, count in row_gene_count.items() if count < 3]
    merged_df.drop(rows_to_drop, axis=0, inplace=True)

    # Save the merged dataframe
    merged_df.to_csv(output_file)
    print("Merged data saved to", output_file)
else:
    print("No dataframes to merge.")
