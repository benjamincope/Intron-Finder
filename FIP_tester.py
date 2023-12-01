import final_intron_finder as FIP


def test_process_fasta():
    sequences = FIP.process_fasta("dummy_FASTA.txt")
    
    for index in range(len(sequences)):
        ID = list(sequences.keys())[index]
        sequence = list(sequences.values())[index]
        
        assert (ID == f"CHR{index + 1}")
        assert (not("\n" in sequence))
    
    
    

def test_process_gff():
    genes = FIP.process_gff("dummy_GFF.txt")
    
    exons = genes[0]
    coding_seqs = genes[1]
    
    for index in range(len(exons)):
        # Everything after the first item is sets of exon-coordinates
        exon_coords = exons[index][1:]
        
        assert exon_coords == coding_seqs[index]
    
    
    
    
def test_check_cds_size():
    genes = FIP.process_gff("dummy_GFF.txt")
    coding_seqs = genes[1]
    print(coding_seqs)
    
    gene_1_coding_seqs = coding_seqs[0] # Sum of CDS's is divisble by 3
    gene_2_coding_seqs = coding_seqs[1] # Sum of CDS's isn't divisible by 3
    
    assert "Success" in FIP.check_cds_size(gene_1_coding_seqs)
    assert "ERROR" in FIP.check_cds_size(gene_2_coding_seqs)
    assert len(coding_seqs) == 6




def test_output_exons():
    sequences = FIP.process_fasta("dummy_FASTA.txt")
    genes = FIP.process_gff("dummy_GFF.txt")
    
    gene_1 = genes[0][0]
    FIP.output_exons(gene_1, sequences)





    
#test_process_fasta()
#test_process_gff()
test_check_cds_size()
#test_output_exons()

