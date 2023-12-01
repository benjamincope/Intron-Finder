## Intron Finder Program - Ben Cope (bec49)
## ----------------------------------------
## To allow graphs to be displayed:
## COMMENT line 178
## UNCOMMENT lines 165-168 and 314-317
## Undo these changes to produce the output file in FASTA format
# -----------------------------------------


import re
import glob
import gzip
import matplotlib.pyplot as plt

home_path = "C:/Users/benja/OneDrive - Aberystwyth University/Desktop/Major Project/"



# Function to return an array of every FASTA and GFF file
def get_file_paths():
    # Locations of directories holding the FASTA and GFF files
    fastas_directory = home_path + "ensemble_fungi/ensemble_fungi/ensembl_fungi_genomes/"
    gffs_directory = home_path + "ensemble_fungi/ensemble_fungi/ensembl_fungi_gff/"
    
    fasta_paths, gff_paths = [], []
    
    # Find every file in the FASTA files directory and the GFF files directory.
    # Files must contain '.fa' / '.gff3' and then optionally '.gz' afterward.
    for path in glob.glob(fastas_directory +"*.fa*"):        
        fasta_paths.append(path)
            
    for path in glob.glob(gffs_directory + "*.gff3*"):
        gff_paths.append(path)
    
    return fasta_paths, gff_paths




def process_fasta(file_name):
    # If the current file is gunzipped, then unzip it
    if (file_name.endswith(".gz")):
        file = gzip.open(file_name, "rt")
    else:
        file = open(file_name, "rt")
        
    data = file.read().split('>')
    file.close()
    del data[0] # Delete empty string
    
    chromosomes = {}
    
    for chrom in data:        
        # Chromosome's ID is the first item on the first line
        chrom_id = chrom.split(" ")[0]
        
        # Chromsome's sequence is everything after the first line
        chrom_sequence = chrom.partition("\n")[2].replace('\n', '')
        
        chromosomes[chrom_id] = chrom_sequence

    return chromosomes




def process_gff(file_name):
    if (file_name.endswith(".gz")):
        file = gzip.open(file_name, "rt")
    else:
        file = open(file_name, "rt")
        
    data = file.read().split('###')
    
    file.close()
    del data[0] # Information before first '###' doesn't represent a DNA sequence
    del data[-1] # Delete empty string at end of file
    
    genome_exons = []
    genome_coding_seqs = []
    row_length = 9 # Number of columns in GFF format
    
    for gene in data:
        # Split the GFF gene on newline and tab characters
        gene = re.split('\n|\t', gene)
        del gene[0], gene[-1] # Remove the empty strings at the start and end
        
        # Will store all exons. First item is information about the gene.
        exons = [gene[0:row_length]]
        
        # Will store coordinates for every CDS
        coding_seqs = []
        
        # Loop through each row of the GFF gene, starting at the 2nd row
        for index in range(row_length, len(gene), row_length):
            if (gene[index + 2] == "exon"):
                # If current row is an exon, add its coordinates to 'exons' array
                exon_coords = [int(gene[index + 3]), int(gene[index + 4])]
                exons.append(exon_coords)
            elif (gene[index + 2] == "CDS"):
                # If current row is a CDS, add its coordinates to 'coding_seqs' array
                cds_coords = [int(gene[index + 3]), int(gene[index + 4])]
                coding_seqs.append(cds_coords)
                
        genome_exons.append(exons)
        genome_coding_seqs.append(coding_seqs)

    return genome_exons, genome_coding_seqs




def check_cds_size(gene_coding_seqs):    
    # Counter to add up size of every CDS
    cds_size = 0
    
    # Iterates through every set of CDS-coordinates for the current gene
    for CDS in gene_coding_seqs:
        # Add size of current CDS, inclusive of its start and stop positions
        cds_size += (CDS[1] - CDS[0]) + 1
        
    # Default outcome message if an error is not caught
    outcome_message = "\n\n***Assembly Success! - CDS of gene is correct\n"
    
    # If coding sequence isn't divisble by 3, there was an issue with splicing
    if ((cds_size % 3) != 0):
        outcome_message = "\n\n***ASSEMBLY ERROR*** - CDS of gene is incomplete\n"
    
    return outcome_message
    

    

def find_introns(gff_gene, fasta_chromosomes):
    # First item of array is an information-array. All other items are exon coords.
    gene_info = gff_gene[0]
    
    # Retrieve all relevant information about the gene
    chromosome_id = gene_info[0]
    attributes = gene_info[8].split(';')
    gene_id = attributes[0].split(':')[1]
    
    # If biotype is tRNA (not protein-coding), then return nothing
    if ("tRNA" in attributes):
        return ""
    
    # Set-up string to store each intron in FASTA format
    intron_info = ""
    
    # Set-up string to store all intron sequences combined together
    combined_intron_seq = ""
    
    # Loop through gene's exons in order to find locations of all introns
    for index in range(2, len(gff_gene)):
        # Intron is from previous exon's end coord to current exon's start coord
        prev_end_coord = gff_gene[index - 1][1]
        current_start_coord = gff_gene[index][0]
        
        intron_coords = [prev_end_coord, current_start_coord]
        size = (intron_coords[1] - intron_coords[0])
        
        seq = fasta_chromosomes[chromosome_id][intron_coords[0] : intron_coords[1]]
        combined_intron_seq += seq
        
        '''
        intron_info += "\n*INTRON*\tCoordinates: " + str(intron_coords) + "\
    \t\tSize: " + str(size) + "\t\tSequence: " + seq
        '''

    
    # Create a new sequence that has a newline at every 60th position, as long
    # as the combined sequence was itleast 200 nucleobases long
    if (len(combined_intron_seq) > 200):
        final_sequence = ""
        for j in range(0, len(combined_intron_seq), 60):
            final_sequence += combined_intron_seq[j : (j + 60)] + "\n"
            
        intron_info += f">{gene_id}_introns\n{final_sequence}"
      
        
        
    # If combined sequence was less than 200 bases long, intron_info will be empty
     
    return intron_info




def output_exons(gff_gene, fasta_chromosomes):
    
    chrom_id = gff_gene[0][0] # ID of chromosome so know where to look in FASTA file
    
    # gff_gene[0] is gene-info so start looping after that item
    for exon_coords in gff_gene[1:]:
        seq = fasta_chromosomes[chrom_id][exon_coords[0] : exon_coords[1]]
        print("\n\n\nCoords:", exon_coords, "\nSequence:", seq)




def display_graphical_statistics(genome_introns):
    # String to be used to combine all intron sequences together
    genome_intron_sequence = ""
    genome_intron_sizes = {}

    lines = genome_introns.split("\n")
    for i in range(len(lines)):
        # Each line starting with *INTRON* has an intron sequence at the end
        if (lines[i].startswith("*INTRON*")):
            intron_sequence = lines[i].split("Sequence: ")[1]
            genome_intron_sequence += intron_sequence
            
            intron_size = int(lines[i].split("Size: ")[1].split("\t")[0])
            if (intron_size in genome_intron_sizes):
                genome_intron_sizes[intron_size] += 1
            else:
                genome_intron_sizes[intron_size] = 1


    
    # ** PLOT 1 **
    # ----------------
    x_markers = ['A', 'C', 'G', 'T']
    y_values = [genome_intron_sequence.count(nt) for nt in x_markers]
    
    
    #plt.subplots_adjust(hspace=1)
    
    plt.subplot(1,2,1)
    plt.bar(x_markers, y_values, width=0.5, color="green")
    plt.title("Nucleotide counts")
    plt.xlabel("Nucleotide")
    plt.ylabel("Occurences in genome (millions)")
    
    
    # ** PLOT 2 **
    # ----------------    
    plt.subplot(1,2,2)
    plt.bar(genome_intron_sizes.keys(), genome_intron_sizes.values())
    plt.xlim([0, 100])
    
    plt.title("Intron sizes")
    plt.xlabel("Size")
    plt.ylabel("Number of introns")
    
    
    plt.subplots_adjust(wspace=0.5) # Stops them overlapping
    plt.show()
    



def main_program():
    # Clear the current introns text file
    intron_file = open("introns_final.txt", 'w')
    intron_file.close()
    
    # Get a list of the paths of every FASTA and GFF file
    fasta_paths, gff_paths = get_file_paths()
    
    # String to hold every genome's introns combined, for purpose of graphs
    all_genome_introns = ""
    
    # Index = how many genomes to run this for
    for index in range(1):
        # Read in the contents of the genome's FASTA and GFF files.
        fasta_chromosomes = process_fasta(fasta_paths[index])
        gff_genes = process_gff(gff_paths[index])
        
        # Find the genome name so can write it to the introns file
        genome_name = fasta_paths[index].split("\\")[1].split(".cds")[0]
        
        # Will store information about all of the introns for this genome
        genome_introns = ""  
        
        # gff_genes[0] is a list of every gene's information + exons
        num_of_genes = len(gff_genes[0])
        
        # Loop through every gene in this genome
        for j in range(num_of_genes):
            # [0] stores the arrays of exons and [1] stores the arrays of CDS's
            gene_exons = gff_genes[0][j]
            gene_coding_seqs = gff_genes[1][j]
            
            # gene_exons stores gene-info, then coords for all exons. So if its length is 2,
            # then there's 1 exon (so no introns) but anything higher means itleast 1 intron
            if (len(gene_exons) > 2):
                # Retrieve information about all of the introns in this gene
                intron_info = find_introns(gene_exons, fasta_chromosomes)
                genome_introns += intron_info
                
            '''
            # Find out if sum of the CDS's for this gene is divisble by 3
            genome_introns += check_cds_size(gene_coding_seqs)
            '''
            
            '''
            # Just run test for the first gene
            if (j == 0):
                output_exons(gene_exons, fasta_chromosomes)
            ''' 

            
        intron_file = open("introns_final.txt", 'a')
        intron_file.write(f"\n\n\n\n{genome_name}\n\n")
        intron_file.write(genome_introns)
        intron_file.close()
        
        
        all_genome_introns += genome_introns
    
        print(f"{index + 1} genomes processed")
    
    '''
    # Output charts to display interesting information about the introns
    display_graphical_statistics(all_genome_introns)    
    ''' 

main_program()






    