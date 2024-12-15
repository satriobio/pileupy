import pysam

def load_refgene(refgene_path):
    """Load reference gene data from file
    
    Args:
        refgene_path: Path to refGene file
        
    Returns:
        List of dictionaries containing gene information
    """
    genes = []
    try:
        with open(refgene_path, 'r') as f:
            for line in f:
                fields = line.strip().split('\t')
                genes.append({
                    'chrom': fields[2],
                    'start': int(fields[4]),
                    'end': int(fields[5]),
                    'name': fields[12],
                    'strand': fields[3],
                    'exon_starts': [int(x) for x in fields[9].strip(',').split(',')],
                    'exon_ends': [int(x) for x in fields[10].strip(',').split(',')]
                })
    except Exception as e:
        print(f"Error reading refGene file: {str(e)}")
        return []
    return genes

def load_bed(bed_path):
    """Load annotation data from BED file
    
    Args:
        bed_path: Path to BED file
        
    Returns:
        List of dictionaries containing annotation information
    """
    annotations = []
    try:
        with open(bed_path, 'r') as f:
            for line in f:
                if line.strip() and not line.startswith('#'):
                    fields = [x for x in line.strip().split() if x]
                    if len(fields) >= 6:
                        annotations.append({
                            'chrom': fields[0],
                            'start': int(fields[1]),
                            'end': int(fields[2]),
                            'name': fields[3],
                            'score': float(fields[4]),
                            'strand': fields[5]
                        })
    except Exception as e:
        print(f"Error reading BED file: {str(e)}")
        return []
    return annotations

def load_cytoband(cytoband_path):
    """Load cytoband data from file
    
    Args:
        cytoband_path: Path to cytoband file
        
    Returns:
        Tuple of (cytobands list, chromosome lengths dict)
    """
    cytobands = []
    chr_lengths = {}
    try:
        with open(cytoband_path, 'r') as f:
            for line in f:
                fields = line.strip().split('\t')
                if len(fields) >= 5:
                    chrom = fields[0]
                    start = int(fields[1])
                    end = int(fields[2])
                    
                    cytobands.append({
                        'chrom': chrom,
                        'start': start,
                        'end': end,
                        'name': fields[3],
                        'stain': fields[4]
                    })
                    
                    chr_lengths[chrom] = max(chr_lengths.get(chrom, 0), end)
    except Exception as e:
        print(f"Error reading cytoband file: {str(e)}")
        return [], {}
    return cytobands, chr_lengths 