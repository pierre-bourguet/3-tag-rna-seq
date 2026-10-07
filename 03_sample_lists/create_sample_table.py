import sys
import glob
import os
from collections import defaultdict
import matplotlib.pyplot as plt

def load_genotype_to_well(table_file):
    genotype_to_well = {}
    with open(table_file, 'r') as file:
        for line in file:
            parts = line.strip().split()
            if len(parts) == 2:
                genotype_to_well[parts[1]] = parts[0]  # Map well to genotype
    return genotype_to_well

def main(directory, table_file, experiment_name):
    # Load mapping of well to genotype
    genotype_to_well = load_genotype_to_well(table_file)
    
    # Reverse the mapping to get genotype from well
    well_to_genotype = {well: genotype for well, genotype in genotype_to_well.items()}
    
    # Find all fastq.gz files
    path_pattern = f"{directory}/*.fastq.gz"
    filenames = glob.glob(path_pattern)
    
    # Replicate tracking and prepare the output
    replicate_count = defaultdict(int)
    sample_output = []
    control_output = [("file", "sample", "size_in_Mb")]
    
    for filename in filenames:
        absolute_path = os.path.abspath(filename)
        basename = os.path.basename(filename)
        well_position = basename.split('_')[0]
        
        # If filename starts with "no_barcode_match", keep the name as it is
        if basename.startswith("no_barcode_match"):
            sample_name = basename.replace('.fastq.gz', '')
        else:
            genotype_description = well_to_genotype.get(well_position, "unknown")
            # Increment and retrieve replicate count
            replicate_count[genotype_description] += 1
            replicate_suffix = f"_R{replicate_count[genotype_description]}"
            sample_name = genotype_description + replicate_suffix if genotype_description != "unknown" else basename.replace('.fastq.gz', '')
        
        file_size_mb = os.path.getsize(filename) / 1024 / 1024  # Convert size to megabytes
        control_output.append((absolute_path, sample_name, file_size_mb))
        
        # If genotype is not unknown and filename doesn't start with "no_barcode_match", add to sample list
        if genotype_description != "unknown" and not basename.startswith("no_barcode_match"):
            sample_output.append(f"{absolute_path}\t{sample_name}\n")

    # Sort control output by file size in descending order while keeping the header at the top
    header = control_output.pop(0)
    control_output.sort(key=lambda x: x[2], reverse=True)
    control_output.insert(0, header)
    
    # Write the results to files
    sample_list_filename = f"sample_list_{experiment_name}.tsv"
    control_file_size_filename = f"control_file_size_{experiment_name}.tsv"

    with open(sample_list_filename, 'w') as file:
        file.writelines(sample_output)
    
    with open(control_file_size_filename, 'w') as file:
        file.write("\t".join(control_output[0]) + "\n")
        for line in control_output[1:]:
            file.write("\t".join(map(str, line)) + "\n")
    
    # Plotting
    sample_names = [item[1] for item in control_output[1:]]  # Skip the header
    file_sizes = [item[2] for item in control_output[1:]]  # Skip the header
    
    fig, ax = plt.subplots(figsize=(10, 20))  # Adjust figure size for better fit and readability
    bars = ax.barh(sample_names, file_sizes, color='skyblue')  # Use horizontal bars
    ax.set_xlabel('File Size (Mb)')
    ax.set_title('Sizes of Demultiplexed Files')
    ax.set_yticks(range(len(sample_names)))
    ax.set_yticklabels(sample_names, rotation=0, ha='right')  # Ensure sample names are readable and right-aligned
    ax.invert_yaxis()  # Largest bar on top

    # Add grid lines at x-axis ticks (every 1000 Mb)
    ax.grid(True, which='major', axis='x', linestyle='--', linewidth=0.5)
    ax.set_xticks(range(0, int(max(file_sizes) + 1000), 1000))  # Set x-ticks every 1000 Mb

    # Adjust the color of y-axis labels based on their content
    y_labels = ax.get_yticklabels()
    for label in y_labels:
        # Use a stricter condition for matching the 'unknown' samples
        if 'unknown' in label.get_text().lower():
            label.set_color('green')  # Set the color of labels for unknown samples to green

    # Remove the white space before the first bar and after the last one
    ax.margins(y=0)  # Set y-axis margins to zero

    plt.tight_layout()  # Automatically adjust subplot parameters to give specified padding
    plt.savefig(f"file_sizes_{experiment_name}.pdf")

    print(f"Generated barplot: file_sizes_{experiment_name}.pdf")

    print(f"Created sample list file: {sample_list_filename}")
    print(f"Created control files list with sizes: {control_file_size_filename}")

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python script.py <directory> <table_file> <experiment_name>")
    else:
        directory, table_file, experiment_name = sys.argv[1], sys.argv[2], sys.argv[3]
        main(directory, table_file, experiment_name)
