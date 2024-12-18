import pysam
import bokeh.plotting as bplt
from bokeh.layouts import column, row
from bokeh.models import ColumnDataSource, Range1d, LinearAxis, Grid, Label, HoverTool, TapTool
from bokeh.models.glyphs import Rect, Text
import numpy as np
from bokeh.io import curdoc
from bokeh.server.server import Server
from bokeh.application.handlers import FunctionHandler
from bokeh.application import Application

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

class Pileupy:
    def __init__(self, region, genome=None, reference=None):
        """Initialize the pileup viewer with basic genome information
        
        Args:
            region: Region string (e.g. 'chr22:24376166-24376456')
            genome: Genome build (default: None) - use either genome or reference
            reference: Path to custom reference file (default: None)
        """
        if genome is None and reference is None:
            raise ValueError("Must provide either genome or reference")
        if genome is not None and reference is not None:
            raise ValueError("Cannot provide both genome and reference")

        self.chrom, pos = region.split(':')
        self.start, self.end = map(int, pos.replace(',','').split('-'))
        
        # Set up genome-specific paths
        self.genome_paths = {
            'hg19': {
                'ref': 'data/hg19.fa',
                'refgene': 'data/ncbiRefSeqSelect.txt',
                'cytoband': 'data/cytoBandIdeo.txt'
            },
            'hg38': {
                'ref': 'hg38.fa',
                'refgene': 'ncbiRefSeqSelect_hg38.txt',
                'cytoband': 'cytoBandIdeo_hg38.txt'
            }
        }
        
        # Initialize core components based on input type
        if genome:
            self.genome = genome
            self.ref = pysam.FastaFile(self.genome_paths[genome]['ref'])
            self.ref_seq = self.ref.fetch(self.chrom, self.start, self.end)
            self.genes = load_refgene(self.genome_paths[genome]['refgene'])
            self.cytobands, self.chr_lengths = load_cytoband(self.genome_paths[self.genome]['cytoband'])
        else:
            self.genome = None
            self.ref = pysam.FastaFile(reference)
            self.ref_seq = self.ref.fetch(self.chrom, self.start, self.end)
            self.genes = []
            self.cytobands = []
            self.chr_lengths = {}
        
        # Initialize shared x range
        self.shared_x_range = Range1d(self.start, self.end)
        
        # Initialize plot components
        self.plots = []
        self.layout = None
        
        # Set up basic styling
        self.read_height = 8
        self.spacing = 2
        self.base_colors = {
            'A': '#00FF00',  # Green
            'T': '#FF0000',  # Red
            'C': '#0000FF',  # Blue
            'G': '#FFA500'   # Orange
        }

        self.filtered_group = []
        
        # Create initial tracks based on input type
        if genome:
            self._create_base_tracks()
        else:
            self.track_reference()
            self._update_layout()

    def _calculate_optimal_height(self, num_reads, min_height=200, max_height=2000):
            """Calculate optimal plot height based on number of reads"""
            reads_per_row = 50  # Match IGV's default
            estimated_rows = num_reads / reads_per_row
            height_per_row = self.read_height + self.spacing  # Match IGV's spacing
            
            calculated_height = max(min_height, min(max_height, estimated_rows * height_per_row))
            return int(calculated_height)

    def on_read_click(self, attr, old, new, source):
        """Callback for read click events"""
        if new:  # If there are selected indices
            reads = []
            for x in new:
                reads.append(source.data['read_names'][x])
                print(f"new {source.data['read_names'][x]}")
            self.filtered_group = reads
        
        # Find and remove existing continuous_multi plots
        # self.plots = [p for p in self.plots if not hasattr(p, 'is_continuous_multi')]
        # self.plots = 
        # # Recreate continuous_multi tracks with updated filtering
        # if hasattr(self, '_continuous_multi_data'):
        #     for data, mode in self._continuous_multi_data:
        #         # Call track_continuous_multi with update_layout=False
        #         self.track_continuous_multi(data, mode, update_layout=False)
        
        # # Update layout once at the end
        # self._update_layout()

    def _create_base_tracks(self):
        """Create the basic genome browser tracks"""
        # Create idiogram track
        self.track_idiogram()
        self.track_reference()
        self.track_genes()
        self._update_layout()

    def track_idiogram(self):
        """Create chromosome idiogram track with cytobands"""
        # Create the figure first with shared x range
        p = bplt.figure(width=800, height=50, tools='', x_range=self.shared_x_range)
        
        # Filter cytobands for current chromosome
        chr_cytobands = [c for c in self.cytobands if c['chrom'] == self.chrom]
        chr_length = self.chr_lengths.get(self.chrom, 0)
        
        if not chr_cytobands or chr_length == 0:
            print("No cytoband data found for chromosome", self.chrom)
            return p
        
        stain_colors = {
            'gneg': '#ffffff',
            'gpos25': '#c0c0c0',
            'gpos50': '#808080',
            'gpos75': '#404040',
            'gpos100': '#000000',
            'acen': '#963232',
            'gvar': '#000000',
            'stalk': '#963232',
        }
        
        # Set up the plot
        p.height = 50
        p.y_range.start = 0
        p.y_range.end = 30
        
        # Fix x_range to show full chromosome, regardless of current view
        p.x_range = Range1d(0, chr_length)
        
        # Draw chromosome backbone
        p.rect(x=chr_length/2,
               y=15,
               width=chr_length,
               height=10,
               fill_color='white',
               line_color='black')
        
        # Draw cytobands
        for band in chr_cytobands:
            width = band['end'] - band['start']
            x = band['start'] + width/2
            
            if band['stain'] == 'acen':
                # Draw centromere as triangle
                if band['name'].startswith('p'):
                    xs = [band['start'], band['end'], band['start']]
                    ys = [10, 15, 20]
                else:
                    xs = [band['start'], band['end'], band['end']]
                    ys = [15, 20, 10]
                
                p.patches(
                    xs=[xs],
                    ys=[ys],
                    fill_color=stain_colors['acen'],
                    line_color='black'
                )
            else:
                # Draw regular band
                p.rect(
                    x=x,
                    y=15,
                    width=width,
                    height=10,
                    fill_color=stain_colors.get(band['stain'], '#ffffff'),
                    line_color='black',
                    line_width=0.5
                )
        
        # Add current view indicator (red line)
        view_pos = (self.start + self.end) / 2
        p.line(
            x=[view_pos, view_pos],
            y=[5, 25],
            line_color='red',
            line_width=2
        )
        
        # Configure axes
        p.xaxis.visible = False
        p.yaxis.visible = False
        p.grid.visible = False
        
        # Disable all interactions
        p.toolbar_location = None
        p.tools = []
        p.toolbar.tools = []
        
        self.plots.append(p)
        self._update_layout()
        return p
    
    def track_reference(self):
        """Create the reference sequence track"""
        # Create the figure first with shared x range
        p = bplt.figure(width=800, height=50, tools='', x_range=self.shared_x_range)
        
        x = np.arange(self.start, self.end)
        y = [2] * len(x)
        
        source = ColumnDataSource({
            'x': x,
            'y': y,
            'text': list(self.ref_seq.upper()),
            'colors': [self.base_colors.get(base.upper(), 'gray') for base in self.ref_seq]
        })
        
        # If zoom level is too far out, show blocks instead of letters
        if (self.end - self.start) > 50:
            p.rect(x='x', y='y', width=1, height=2,
                   fill_color='colors', line_color=None, source=source)
        else:
            p.text(x='x', y='y', text='text', text_color='colors', source=source,
                   text_baseline="middle", text_align="center")
        
        # Configure axes
        p.xaxis.visible = True
        p.xaxis.axis_label_text_font_size = '10pt'
        p.xaxis[0].formatter.use_scientific = False
        p.xaxis.axis_line_color = 'black'
        p.xaxis.major_tick_line_color = 'black'
        p.xaxis.minor_tick_line_color = None
        p.above = p.below  # Move axis to top
        p.below = []  # Remove bottom axis
        
        # Remove y-axis completely
        p.yaxis.visible = False

        # Remove all grid lines
        p.grid.visible = False
        
        
        self.plots.append(p)
        self._update_layout()
        return p

    def track_genes(self):
        """Create gene track visualization"""
        # Create the figure first with shared x range
        p = bplt.figure(width=800, height=50, tools='', x_range=self.shared_x_range)
        
        # Filter genes in view range
        visible_genes = [g for g in self.genes 
                        if g['chrom'] == self.chrom and 
                        not (g['end'] < self.start or g['start'] > self.end)]

        if not visible_genes:
            # Add placeholder text if no genes are found
            p.text(x=[self.start + (self.end - self.start)/2],
                   y=[25],
                   text=['No genes in this region'],
                   text_font_size='8pt',
                   text_baseline="middle",
                   text_align="center",
                   text_color='gray')
        else:
            # Draw genes
            for gene in visible_genes:
                y_pos = 25  # Base y position
                exon_height = 1.5  # Reduced exon height
                line_height = 1  # Height of the connecting line
                arrow_size = 1  # Size of direction arrows

                # Draw gene body (thin blue line)
                p.line(x=[gene['start'], gene['end']],
                       y=[y_pos, y_pos],
                       line_color='blue',
                       line_width=line_height)

                # Draw direction arrows along the gene body
                arrow_spacing = 50  # Pixels between arrows
                gene_length = gene['end'] - gene['start']
                num_arrows = max(1, min(10, gene_length // arrow_spacing))
                arrow_width = 6  # Width of the arrow
                arrow_height = 1  # Height of the arrow

                for i in range(num_arrows):
                    arrow_x = gene['start'] + (i + 0.5) * (gene_length / num_arrows)
                    if gene['strand'] == '+':
                        # Right-pointing arrow
                        xs = [
                            arrow_x - arrow_width/2,  # Left point
                            arrow_x + arrow_width/2,  # Right point (tip)
                            arrow_x - arrow_width/2   # Back to start
                        ]
                        ys = [
                            y_pos - arrow_height/2,  # Bottom point
                            y_pos,                   # Middle point (tip)
                            y_pos + arrow_height/2   # Top point
                        ]
                    else:  # '-' strand
                        # Left-pointing arrow
                        xs = [
                            arrow_x + arrow_width/2,  # Right point
                            arrow_x - arrow_width/2,  # Left point (tip)
                            arrow_x + arrow_width/2   # Back to start
                        ]
                        ys = [
                            y_pos - arrow_height/2,  # Bottom point
                            y_pos,                   # Middle point (tip)
                            y_pos + arrow_height/2   # Top point
                        ]
                    
                    p.patches(
                        xs=[xs],
                        ys=[ys],
                        fill_color='white',
                        line_color='black'
                    )

                # Draw exons (blue rectangles)
                for start, end in zip(gene['exon_starts'], gene['exon_ends']):
                    if not (end < self.start or start > self.end):  # Only draw visible exons
                        p.rect(x=[(start + end)/2],
                               y=[y_pos],
                               width=[end - start],
                               height=[exon_height],
                               color='blue',
                               line_color='black',
                               alpha=0.8)

                # Add gene name in the middle of the gene, level with exons
                p.text(x=[(gene['start'] + gene['end']) / 2],
                       y=[y_pos],  # Position text at the same level as exons
                       text=[gene['name']],
                       text_font_size='8pt',
                       text_baseline="middle",
                       text_align="center")

        # Configure axes
        p.xaxis.visible = False
        p.yaxis.visible = False
        p.grid.visible = False

        p.toolbar_location = None

        self.plots.append(p)
        self._update_layout()
        return p

    def track_coverage(self, bam):
        """Create the coverage track"""
        # Create the figure first with shared x range
        p = bplt.figure(width=800, height=50, tools='', x_range=self.shared_x_range)
        
        coverage = bam.count_coverage(self.chrom, self.start, self.end)
        total_coverage = np.sum(coverage, axis=0)
        
        source = ColumnDataSource({
            'x': np.arange(self.start, self.end),
            'y': total_coverage,
            'height': total_coverage
        })
        
        p.vbar(x='x', top='y', width=1, source=source,
               fill_color='gray', line_color='gray', line_width=5)
        
        p.xaxis.visible = False
        min_val = 0
        print(total_coverage)
        max_val = int(max(total_coverage))
        p.yaxis.ticker = [min_val, max_val]

        p.grid.visible = False

        p.toolbar_location = None
        
        self.plots.append(p)
        self._update_layout()
        return p

    def track_alignment(self, bam_path, mode='collapsed'):
        """Add an alignment track to the visualization
        
        Args:
            bam_path: Path to BAM file
            mode: Display mode ('collapsed' or 'expanded')
        """
        bam = pysam.AlignmentFile(bam_path, "rb")
        reads = list(bam.fetch(self.chrom, self.start, self.end))
        
        # Create new figure with shared x range
        p = bplt.figure(width=800, height=200, 
                        tools='xwheel_zoom,xwheel_pan',
                        x_range=self.shared_x_range)
        
        # Calculate optimal height ONCE using the correct method name
        optimal_height = self._calculate_optimal_height(len(reads))
        p.height = optimal_height
        p.y_range = Range1d(0, optimal_height)
        
        # Draw alignment track
        print(f"Number of reads found: {len(reads)}")
        
        read_data = {
            'xs': [], 'ys': [], 'colors': [],
            'snp_xs': [], 'snp_ys': [], 'snp_colors': [],  # Added snp_colors
            'snp_heights': []  # Added for line heights
        }
        
        region_length = self.end - self.start
        max_rows = optimal_height // (self.read_height + self.spacing)
        occupied = np.zeros((max_rows, region_length), dtype=bool)
        
        processed_count = 0
        for read_idx, read in enumerate(reads):
            try:
                if not read.is_unmapped and read.reference_length:
                    read_start_rel = read.reference_start - self.start
                    read_end_rel = read_start_rel + read.reference_length
                    
                    # Find available row
                    row = 0
                    while row < max_rows:
                        start_idx = max(0, read_start_rel)
                        end_idx = min(region_length, read_end_rel)
                        check_start = max(0, start_idx - 1)
                        check_end = min(region_length, end_idx + 1)
                        
                        if start_idx < end_idx and not any(occupied[row][check_start:check_end]):
                            occupied[row][start_idx:end_idx] = True
                            break
                        row += 1
                    
                    if row < max_rows:
                        processed_count += 1
                        y_pos = optimal_height - (row + 1) * (self.read_height + self.spacing)
                        
                        # Calculate arrow points
                        x_start = read.reference_start
                        x_end = read.reference_start + read.reference_length
                        y_top = y_pos + self.read_height
                        y_bottom = y_pos
                        arrow_width = 2  # Changed from 3 to 2
                        
                        y_mid = y_pos + (self.read_height / 2)  # Calculate middle y-position
                        
                        if read.is_reverse:
                            # Left-pointing arrow (reverse strand)
                            xs = [
                                x_start, x_end,          # Bottom line
                                x_end, x_end,            # Right edge
                                x_end, x_start + arrow_width,  # Top line
                                x_start + arrow_width, x_start + arrow_width/2,  # Arrow point top (less sharp)
                                x_start + arrow_width/2, x_start + arrow_width,  # Arrow point bottom (less sharp)
                                x_start + arrow_width              # Back to start
                            ]
                            ys = [
                                y_bottom, y_bottom,      # Bottom line
                                y_bottom, y_top,         # Right edge
                                y_top, y_top,           # Top line
                                y_top, y_mid,           # Arrow point top
                                y_mid, y_bottom,        # Arrow point bottom
                                y_bottom                # Back to start
                            ]
                        else:
                            # Right-pointing arrow (forward strand)
                            xs = [
                                x_start, x_end - arrow_width,  # Bottom line
                                x_end - arrow_width, x_end - arrow_width/2,    # Arrow point bottom (less sharp)
                                x_end - arrow_width/2, x_end - arrow_width,    # Arrow point top (less sharp)
                                x_end - arrow_width, x_start,  # Back to start
                                x_start, x_start              # Left edge
                            ]
                            ys = [
                                y_bottom, y_bottom,      # Bottom line
                                y_bottom, y_mid,         # Arrow point
                                y_mid, y_top,           # Top line
                                y_top, y_top,           # Back to start
                                y_top, y_bottom         # Left edge
                            ]
                        
                        read_data['xs'].append(xs)
                        read_data['ys'].append(ys)
                        read_data['colors'].append('#E0E0E0')
                        
                        # Check for SNPs
                        if read.query_sequence and self.ref_seq:
                            ref_pos = read.reference_start - self.start
                            query_pos = 0

                            if read.query_name == 'HWI-BRUNOP16X_0001:7:61:2300:181533#0':
                                print(length)
                                print(read.cigartuples)

                                for op, length in read.cigartuples:
                                    print(op)
                                    print(length)
                                    if op == 0:  # Match/mismatch
                                        for i in range(length):

                                            print(ref_pos, i)
                                            print(ref_pos, self.ref_seq)
                                            print(query_pos, len(read.query_sequence))

                                            if (ref_pos + i >= 0 and 
                                                ref_pos + i < len(self.ref_seq) and 
                                                query_pos + i < len(read.query_sequence)):
                                                
                                                ref_base = self.ref_seq[ref_pos + i].upper()
                                                query_base = read.query_sequence[query_pos + i].upper()
                                                
                                                print('here')

                                                if ref_base != query_base:
                                                    # print('HERE')
                                                    read_data['snp_xs'].append(read.reference_start + i)
                                                    read_data['snp_ys'].append(y_pos)
                                                    read_data['snp_colors'].append('#FF00FF')
                                                    read_data['snp_heights'].append(self.read_height)
                                        
                                        ref_pos += length
                                        query_pos += length
                                    elif op == 1:  # Insertion
                                        query_pos += length
                                    elif op == 2:  # Deletion
                                        ref_pos += length
                            else:
                            
                                for op, length in read.cigartuples:
                                    if op == 0:  # Match/mismatch
                                        for i in range(length):
                                            if (ref_pos + i >= 0 and 
                                                ref_pos + i < len(self.ref_seq) and 
                                                query_pos + i < len(read.query_sequence)):
                                                
                                                ref_base = self.ref_seq[ref_pos + i].upper()
                                                query_base = read.query_sequence[query_pos + i].upper()
                                                
                                                if ref_base != query_base:
                                                    read_data['snp_xs'].append(read.reference_start + i)
                                                    read_data['snp_ys'].append(y_pos)
                                                    read_data['snp_colors'].append(self.base_colors.get(query_base, '#808080'))
                                                    read_data['snp_heights'].append(self.read_height)
                                        
                                        ref_pos += length
                                        query_pos += length
                                    elif op == 1:  # Insertion
                                        query_pos += length
                                    elif op == 2:  # Deletion
                                        ref_pos += length
                                    
            except Exception as e:
                print(f"Error processing read {read_idx}: {str(e)}")
                continue
        print(f"\nProcessed {processed_count} out of {len(reads)} reads")
        
        if len(read_data['xs']) > 0:
            # Add additional data for tooltips
            read_data['read_names'] = []
            read_data['mapping_qualities'] = []
            read_data['cigar_strings'] = []
            read_data['read_lengths'] = []
            
            # Populate the additional data while processing reads
            for read in reads:
                if not read.is_unmapped and read.reference_length:
                    read_data['read_names'].append(read.query_name)
                    read_data['mapping_qualities'].append(read.mapping_quality)
                    read_data['cigar_strings'].append(read.cigarstring)
                    read_data['read_lengths'].append(len(read.query_sequence) if read.query_sequence else 0)
            
            source = ColumnDataSource(read_data)
            
            # Add tooltips
            tooltips = [
                ('Read Name', '@read_names'),
                ('Mapping Quality', '@mapping_qualities'),
                ('CIGAR', '@cigar_strings'),
                ('Length', '@read_lengths')
            ]
            
            # Add HoverTool
            hover_tool = HoverTool(tooltips=tooltips, renderers=[])
            p.add_tools(hover_tool)
            
            # Draw reads as interactive patches
            patches = p.patches(
                xs='xs',
                ys='ys',
                source=source,
                fill_color='colors',
                line_color=None,
                alpha=0.8,
                name='read_patches'  # Add name for identification
            )
            
            # Add hover tool to the patches renderer
            hover_tool.renderers.append(patches)
            
            # Add tap (click) callback
            source.selected.on_change('indices', lambda attr, old, new: self.on_read_click(attr, old, new, source))
            
            # Add TapTool
            p.add_tools(TapTool())
            
            # Draw SNPs
            if len(read_data['snp_xs']) > 0:
                snp_source = ColumnDataSource({
                    'x': read_data['snp_xs'],
                    'y': read_data['snp_ys'],
                    'y1': [y + h for y, h in zip(read_data['snp_ys'], read_data['snp_heights'])],
                    'color': read_data['snp_colors']
                })
                
                p.segment(
                    x0='x', y0='y',
                    x1='x', y1='y1',  # Use pre-calculated y1 values
                    color='color',
                    line_width=2,
                    source=snp_source,
                    alpha=0.8
                )
        
        # Adjust y_range based on actual data
        if len(read_data['ys']) > 0:
            min_y = min([min(y) for y in read_data['ys']]) - self.read_height
            max_y = optimal_height
            p.y_range = Range1d(min_y, max_y)
        
        # Remove x-axis ticks and labels
        p.xaxis.visible = False
        p.yaxis.visible = False
        p.grid.visible = False
        p.toolbar.logo = None

        p2 = self.track_coverage(bam)
        
        self.plots.append(p2)
        self.plots.append(p)
        self._update_layout()
        return p

    def track_annotation(self, bed_path):
        """Add an annotation track from BED file
        
        Args:
            bed_path: Path to BED file
        """
        # Create new figure
        p = bplt.figure(width=800, height=50, tools='', x_range=self.shared_x_range)
        
        # Load and filter annotations
        annotations = load_bed(bed_path)
        visible_annotations = [f for f in annotations 
                             if f['chrom'] == self.chrom and 
                             not (f['end'] < self.start or f['start'] > self.end)]
        
        if not visible_annotations:
            p.text(x=[self.start + (self.end - self.start)/2],
                  y=[1],
                  text=['No annotations found'],
                  text_font_size='8pt',
                  text_baseline="middle",
                  text_align="center",
                  text_color='gray')
        else:
            # Draw annotations
            for annotation in visible_annotations:
                p.rect(x=[(annotation['start'] + annotation['end'])/2],
                      y=[1],
                      width=[annotation['end'] - annotation['start']],
                      height=[10],
                      color='red',
                      alpha=0.7)
                
                p.text(x=[(annotation['start'] + annotation['end'])/2],
                      y=[1],
                      text=[annotation['name']],
                      text_font_size='8pt',
                      text_baseline="middle",
                      text_align="center")
        
        p.xaxis.visible = False
        p.yaxis.visible = False
        p.grid.visible = False

        p.toolbar_location = None

        self.plots.append(p)
        self._update_layout()
        return p

    def track_continuous(self, df, mode='line'):
        """Add a continuous data track from pandas DataFrame
        
        Args:
            df: Pandas DataFrame with required columns: chrom, start, end
                Optional column: signal (defaults to 1 if not provided)
            mode: Plot type ('line', 'scatter', or 'bar')
        """
        # Create new figure
        p = bplt.figure(width=800, height=70, tools='', x_range=self.shared_x_range)
        
        # Filter data for current chromosome and view range
        mask = (df['chrom'] == self.chrom) & \
               (df['start'] <= self.end) & \
               (df['end'] >= self.start)
        visible_data = df[mask].copy()
        
        if visible_data.empty:
            # Show placeholder text if no data
            p.text(x=[self.start + (self.end - self.start)/2],
                   y=[25],
                   text=['No data in this region'],
                   text_font_size='8pt',
                   text_baseline="middle",
                   text_align="center",
                   text_color='gray')
        else:
            # Use 'signal' column if it exists, otherwise default to 1
            if 'signal' not in visible_data.columns:
                visible_data['signal'] = 1
            
            source = ColumnDataSource({
                'x': visible_data['start'],
                'y': visible_data['signal']
            })
            
            # Plot according to mode
            if mode == 'bar':
                p.vbar(x='x', top='y', bottom=0, width=1, source=source,
                       color='navy', alpha=0.6)
            elif mode == 'line':
                p.line(x='x', y='y', source=source,
                       line_width=2, color='navy', alpha=0.6)
            elif mode == 'scatter':
                p.scatter(x='x', y='y', source=source,
                         size=8, color='navy', alpha=0.6)
            
            # Add hover tool
            hover = HoverTool(tooltips=[
                ('Position', '@x'),
                ('Value', '@y{0.000}')
            ])
            p.add_tools(hover)
        
        # Configure axes
        p.xaxis.visible = False
        p.yaxis.visible = True
        p.grid.visible = False
        p.toolbar.logo = None
        
        self.plots.append(p)
        self._update_layout()
        return p

    def track_continuous_multi(self, df, mode='bar', update_layout=True):
        """Add a continuous data track from pandas DataFrame
        
        Args:
            df: Pandas DataFrame with columns: chrom, start, end, read_name, signal
            mode: Plot type ('bar', 'line', or 'scatter')
            update_layout: Whether to update the layout after adding the track
        """
        # Store the original data and settings for later updates
        if not hasattr(self, '_continuous_multi_data'):
            self._continuous_multi_data = []
        self._continuous_multi_data.append((df.copy(), mode))
        
        # Create new figure
        p = bplt.figure(width=800, height=70, tools='', x_range=self.shared_x_range)
        
        # Mark this plot as a continuous_multi plot
        # p.is_continuous_multi = True
        
        # Filter data for current chromosome and view range
        mask = (df['chrom'] == self.chrom) & \
               (df['start'] <= self.end) & \
               (df['end'] >= self.start)
        visible_data = df[mask].copy()
        
        # Apply read name filtering if self.filtered_group is not empty
        if self.filtered_group:
            visible_data = visible_data[visible_data['read_name'].isin(self.filtered_group)]
        
        if visible_data.empty:
            # Show placeholder text if no data
            p.text(x=[self.start + (self.end - self.start)/2],
                   y=[25],
                   text=['No data in this region'],
                   text_font_size='8pt',
                   text_baseline="middle",
                   text_align="center",
                   text_color='gray')
        else:
            # Get unique read names and assign y-positions
            unique_reads = visible_data['read_name'].unique()
            spacing = 1
            y_positions = {read: (i + 1) * spacing for i, read in enumerate(unique_reads)}
            
            # Add y-position to the data
            visible_data['y_pos'] = visible_data['read_name'].map(y_positions)
            
            # Plot each read separately with its own color
            colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd']  # Default Bokeh colors
            
            for i, read_name in enumerate(unique_reads):
                read_data = visible_data[visible_data['read_name'] == read_name]
                color = colors[i % len(colors)]  # Cycle through colors
                
                source = ColumnDataSource({
                    'x': read_data['start'],
                    'y': read_data['signal'], # 'y': read_data['signal'] + read_data['y_pos'],  # Offset signal by y_position
                    'base_y': read_data['y_pos'],  # Base y position for this read
                    'read_name': read_data['read_name'],
                    'signal': read_data['signal']
                })
                
                # Plot according to mode
                if mode == 'bar':
                    p.vbar(x='x', top='y', bottom='base_y', width=1, source=source,
                           color=color, alpha=0.6)
                elif mode == 'line':
                    p.line(x='x', y='y', source=source,
                           line_width=2, color=color, alpha=0.6)
                elif mode == 'scatter':
                    p.scatter(x='x', y='y', source=source,
                             size=8, color=color, alpha=0.6)
            
            # Add hover tool
            hover = HoverTool(tooltips=[
                ('Read Name', '@read_name'),
                ('Position', '@x'),
                ('Signal', '@signal{0.000}')  # Show original signal value
            ])
            p.add_tools(hover)
            
            # Configure legend
            # p.legend = False
            # p.legend.click_policy = "hide"
            # p.legend.location = "right"
            # p.legend.background_fill_alpha = 0.7
            
            # Set y-range to accommodate all reads plus some padding
            # max_y = max(visible_data['signal'].max() + len(unique_reads) * spacing)
            # p.y_range.start = 0
            # p.y_range.end = max_y + spacing
        
        # Configure axes
        p.xaxis.visible = False
        p.yaxis.visible = True
        p.grid.visible = False
        p.toolbar.logo = None
        
        self.plots.append(p)
        
        # Only update layout if requested
        if update_layout:
            self._update_layout()
        
        return p

    def _update_layout(self):
        """Update the layout after adding new tracks"""
        if not self.plots:
            print("Warning: No plots to create layout")
            return
        
        try:
            self.layout = column(self.plots)
            print(f"Layout updated with {len(self.plots)} plots")
        except Exception as e:
            print(f"Error updating layout: {str(e)}")

    def show(self):
        """Display the visualization"""
        bplt.show(self.layout)

def setup(doc):
    """Setup function for Bokeh server"""
    try:
        print("Starting setup...")
        
        # Initialize browser
        browser = Pileupy('chr22:24376166-24376456', genome='hg19')
        print("Browser initialized")
        
        # Add tracks
        print("Adding tracks...")
        # browser.track_alignment('data/gstt1_sample.bam')
        print("Alignment track added")
        # browser.track_annotation('data/mod.bed')
        print("Annotation track added")
        
        if browser.layout is None:
            print("Warning: Browser layout is None")
            return
            
        print(f"Setup complete. Layout has {len(browser.plots)} plots")
        doc.add_root(browser.layout)
        
    except Exception as e:
        print(f"Error in setup: {str(e)}")
        import traceback
        traceback.print_exc()
        # Show error in the browser
        p = bplt.figure(width=800, height=400)
        p.text(x=0, y=0, text=[f"Error: {str(e)}"])
        doc.add_root(column([p]))

def modify_doc(doc):
    """This function will be called by the bokeh server"""
    # Initialize the viewer with default parameters
    print("Starting setup...")
        
    # Initialize browser
    browser = Pileupy('chr22:24376166-24376456', genome='hg19')
    print("Browser initialized")
    
    # Add tracks
    print("Adding tracks...")
    browser.track_alignment('data/gstt1_sample.bam')
    print("Alignment track added")
    browser.track_annotation('data/mod.bed')
    print("Annotation track added")

    # Create mock continuous data
    import pandas as pd
    import numpy as np
    
    # Define specific read names
    read_names = [
        "HWI-BRUNOP16X_0001:7:27:7557:21734#0",
        "HWI-BRUNOP16X_0001:7:42:7486:110377#0",
        "HWI-BRUNOP16X_0001:7:5:13679:157971#0",
        "HWI-BRUNOP16X_0001:7:61:16093:161756#0",
        "HWI-BRUNOP16X_0001:7:7:13539:38639#0"
    ]
    
    # Generate continuous data for each read
    mock_data = []
    view_range = browser.end - browser.start
    points_per_read = 100  # Increased for smoother signals
    
    for i, read_name in enumerate(read_names):
        # Create evenly spaced points across the view range
        positions = np.linspace(browser.start, browser.end, points_per_read)
        
        # Generate unique random signal for each read
        # Using different parameters for each read to create variety
        mean = np.random.uniform(0.3, 0.7)  # Random mean between 0.3 and 0.7
        std = np.random.uniform(0.1, 0.2)    # Random standard deviation
        
        # Generate base signal
        base_signal = np.random.normal(mean, std, points_per_read)
        
        # Add some smooth variations using sine waves
        x = np.linspace(0, 4*np.pi, points_per_read)
        wave = 0.2 * np.sin(x + i) + 0.1 * np.sin(2*x + i*2)
        
        # Combine base signal with wave and ensure all values are positive
        signal = np.maximum(0, base_signal + wave)
        
        for pos, val in zip(positions, signal):
            mock_data.append({
                'chrom': 'chr22',
                'start': int(pos),
                'end': int(pos) + 10,
                'read_name': read_name,
                'signal': val
            })
    
    # Create DataFrame
    df = pd.DataFrame(mock_data)
    
    # Sort by position and read_name
    df = df.sort_values(['read_name', 'start'])
    
    # Add the continuous tracks with mock data
    # browser.track_continuous_multi(df, mode='bar')
    browser.track_continuous_multi(df, mode='scatter')
    browser.track_continuous_multi(df, mode='line')
    print("Continuous tracks added")
    
    # Generate mock methylation data
    methylation_data = []
    positions = np.linspace(browser.start, browser.end, 50)  # 50 methylation sites
    
    # Create base methylation pattern with some sites highly methylated and others not
    base_pattern = np.concatenate([
        np.random.uniform(0.8, 1.0, 15),    # Highly methylated region
        np.random.uniform(0.1, 0.3, 20),    # Low methylation region
        np.random.uniform(0.7, 0.9, 15)     # Another highly methylated region
    ])
    
    # Add some noise to make it more realistic
    methylation_values = base_pattern + np.random.normal(0, 0.05, len(base_pattern))
    methylation_values = np.clip(methylation_values, 0, 1)  # Ensure values stay between 0 and 1
    
    for pos, meth in zip(positions, methylation_values):
        methylation_data.append({
            'chrom': 'chr22',
            'start': int(pos),
            'end': int(pos) + 10,
            'read_name': 'methylation',
            'signal': meth
        })
    
    # Create DataFrame for methylation data
    meth_df = pd.DataFrame(methylation_data)
    meth_df = meth_df.sort_values('start')
    
    # Add methylation track
    browser.track_continuous(meth_df)
    print("Methylation track added")

    # Add the layout to the document
    doc.add_root(browser.layout)
    doc.title = "Pileup Viewer"

# Remove the __main__ block and replace with server setup
if __name__ == '__main__':
    # Option 1: Run with 'bokeh serve'
    # modify_doc(curdoc())
    
    # Option 2: Run standalone (uncomment below if you want to run without bokeh serve command)
    """
    """
    def bk_worker():
        server = Server({'/': modify_doc}, num_procs=1)
        server.start()
        server.io_loop.add_callback(server.show, "/")
        server.io_loop.start()

    print('Opening Bokeh application on http://localhost:5006/')
    bk_worker()