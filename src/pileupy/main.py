from .data import load_refgene, load_bed, load_cytoband

import pysam
import bokeh.plotting as bplt
from bokeh.layouts import column, row
from bokeh.models import ColumnDataSource, Range1d, LinearAxis, Grid, Label, HoverTool, TapTool, Select, TextInput, Button, Div
from bokeh.models.glyphs import Rect, Text
import numpy as np
from bokeh.io import curdoc
from bokeh.server.server import Server
from bokeh.application.handlers import FunctionHandler
from bokeh.application import Application
import time  # Add this import at the top of the file

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
        else:
            self.reference = reference

        
        # Initialize shared x rangeg
        
        
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

        if self.genome:
            self.ref = pysam.FastaFile(self.genome_paths[self.genome]['ref'])
            self.ref_seq = self.ref.fetch(self.chrom, self.start, self.end)
            self.genes = load_refgene(self.genome_paths[self.genome]['refgene'])
            self.cytobands, self.chr_lengths = load_cytoband(self.genome_paths[self.genome]['cytoband'])
        else:
            self.genome = None
            self.ref = pysam.FastaFile(self.reference)
            self.genes = []
            self.cytobands = []
            self.chr_lengths = {}
        
        self.shared_x_range = Range1d(self.start, self.end)

        self.track_control()
        self.track_idiogram()
        self.track_reference()
        self.track_genes()
        self._update_layout()

    def track_control(self):
        """Create control panel with navigation elements"""
        # Create chromosome dropdown
        available_chroms = sorted(list(self.chr_lengths.keys()), key=_chr_sort_key)

        chrom_select = Select(
            value=self.chrom,
            options=available_chroms,
            width=100
        )

        # Create position input fields
        start_input = TextInput(
            value=str(self.start),
            width=120
        )
        
        end_input = TextInput(
            value=str(self.end),
            width=120
        )

        # Create zoom buttons
        zoom_out = Button(label="-", button_type="primary", width=15)
        zoom_in = Button(label="+", button_type="primary", width=15)
        
        # Create location text display
        # loc_text = Div(
        #     text=f"{self.chrom}:{self.start:,}-{self.end:,}",
        #     width=200
        # )

        def redraw_tracks():
            """Redraw all tracks with updated coordinates"""
            # Store existing addon tracks
            addons = self.plots[4:]
            
            # Clear existing plots and recreate base tracks
            self.plots = []
            self._create_base_tracks()
            
            # Process addon tracks in batch
            if addons:
                for addon in addons:
                    if len(addon.tags) > 0 and addon.tags[0].get('type') == 'alignment':
                        self.track_alignment(addon.tags[0].get('bam_path'))
                    if len(addon.tags) > 0 and addon.tags[0].get('type') == 'annotation':
                        self.track_alignment(addon.tags[0].get('bed_path'))
            self._update_layout()

            if curdoc().session_context:
                curdoc().clear()
                curdoc().add_root(self.layout)

        def update_region():
            """Update the viewed region"""
            try:
                new_start = int(start_input.value.replace(',', ''))
                new_end = int(end_input.value.replace(',', ''))
                    
                # Update coordinates
                self.start = new_start
                self.end = new_end
                
                # Update location text
                # loc_text.text = f"{self.chrom}:{self.start:,}-{self.end:,}"
                
                # Redraw tracks with new chromosome
                redraw_tracks()
                
            except ValueError as e:
                print(f"Invalid position: {str(e)}")

        def update_chromosome(attr, old, new):
            """Handle chromosome selection"""
            self.chrom = new
            # Reset view to full chromosome
            self.start = 0
            self.end = self.chr_lengths.get(new, 0)
            
            # Update inputs and location text
            start_input.value = str(self.start)
            end_input.value = str(self.end)
            # loc_text.text = f"{self.chrom}:{self.start:,}-{self.end:,}"
            
            # Redraw tracks with new chromosome
            redraw_tracks()

        def zoom(direction):
            """Handle zoom in/out"""
            center = (self.start + self.end) // 2
            current_width = self.end - self.start
            
            if direction == 'in':
                new_width = max(20, current_width // 2)  # Don't zoom in closer than 100bp
            else:
                new_width = min(current_width * 2, self.chr_lengths.get(self.chrom, 0))
                
            half_width = new_width // 2
            new_start = max(0, center - half_width)
            new_end = min(self.chr_lengths.get(self.chrom, 0), center + half_width)
            
            # Update values
            self.start = new_start
            self.end = new_end
            
            # Update inputs and location text
            start_input.value = str(new_start)
            end_input.value = str(new_end)
            # loc_text.text = f"{self.chrom}:{new_start:,}-{new_end:,}"
            
            # Redraw tracks with new zoom level
            redraw_tracks()

        # Connect callbacks
        chrom_select.on_change('value', update_chromosome)
        start_input.on_change('value', lambda attr, old, new: update_region())
        end_input.on_change('value', lambda attr, old, new: update_region())
        zoom_in.on_click(lambda: zoom('in'))
        zoom_out.on_click(lambda: zoom('out'))

        # Create control panel layout
        controls = row(
            chrom_select,
            start_input,
            end_input,
            zoom_out,
            zoom_in,
            spacing=10,
            margin=(0, 0, 0, 20)
        )

        controls_layout = column(controls, sizing_mode="stretch_width")
        self.plots.append(controls_layout)        
        self._update_layout()
        
        return controls_layout
    
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
        p.rect(
            x=(self.start + self.end) / 2,  # center x position
            y=15,  # center y position (same as chromosome backbone)
            width=self.end - self.start,  # width spans the current view
            height=10,  # same height as chromosome backbone
            fill_color='red',
            fill_alpha=0.2,
            line_color='red',
            line_width=1
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
         
        # If zoom level is too far out, disable blocks
        if (self.end - self.start) >= 10000:
            p.text(x=[self.start + (self.end - self.start)/2],
                    y=[25],
                    text=['Zoom in to see reference sequence'],
                    text_font_size='8pt',
                    text_baseline="middle",
                    text_align="center",
                    text_color='gray')
        else:
            self.ref_seq = self.ref.fetch(self.chrom, self.start, self.end)
            source = ColumnDataSource({
                'x': x,
                'y': y,
                'text': list(self.ref_seq.upper()),
                'colors': [self.base_colors.get(base.upper(), 'gray') for base in self.ref_seq]
            })
            # If zoom level is too far out, show blocks instead of letters
            if ((self.end - self.start) < 10000) & ((self.end - self.start) >= 50):
                p.rect(x='x', y='y', width=1, height=2,
                    fill_color='colors', line_color=None, source=source)
            # If zoom level is close enough, show letters
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
        p.y_range = Range1d(0, 50)

        # Time the filtering operation
        filter_start = time.time()
        visible_genes = [g for g in self.genes 
                        if g['chrom'] == self.chrom and 
                        not (g['end'] < self.start or g['start'] > self.end)]
        filter_time = time.time() - filter_start

        # Prepare batch data for gene bodies
        gene_lines_x = []
        gene_lines_y = []
        
        # Prepare batch data for exons
        exon_xs = []
        exon_ys = []
        exon_widths = []
        
        # Prepare batch data for gene names
        text_xs = []
        text_ys = []
        text_labels = []

        y_pos = 35
        exon_height = 20
        line_height = 1

        for gene in visible_genes:
            # Add gene body line coordinates
            if (self.end - self.start) < 1000000:
                gene_lines_x.extend([gene['start'], gene['end'], None])  # None creates a break between segments
                gene_lines_y.extend([y_pos, y_pos, None])

            # Add exon coordinates
            for start, end in zip(gene['exon_starts'], gene['exon_ends']):
                if not (end < self.start or start > self.end):
                    exon_xs.append((start + end)/2)
                    exon_ys.append(y_pos)
                    exon_widths.append(end - start)

            if (self.end - self.start) < 1000000:
                text_xs.append((gene['start'] + gene['end']) / 2)
                text_ys.append(15)
                text_labels.append(gene['name'])

        # Batch plot gene body lines
        if gene_lines_x:
            p.line(x=gene_lines_x, y=gene_lines_y, line_color='blue', line_width=line_height)

        # Batch plot exons
        if exon_xs:
            p.rect(x=exon_xs, y=exon_ys, width=exon_widths, height=exon_height,
                   color='blue', line_color=None)

        # Batch plot gene names
        if text_xs:
            p.text(x=text_xs, y=text_ys, text=text_labels,
                   text_font_size='8pt', text_baseline="middle", text_align="center")
   
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
        p.xaxis.visible = False
        p.yaxis.visible = False
        if (self.end - self.start) >= 10000:
            p.text(x=[self.start + (self.end - self.start)/2],
                    y=[25],
                    text=['Zoom in to view coverage'],
                    text_font_size='8pt',
                    text_baseline="middle",
                    text_align="center",
                    text_color='gray')
        else:
            coverage = bam.count_coverage(self.chrom, self.start, self.end)
            total_coverage = np.sum(coverage, axis=0)
            
            source = ColumnDataSource({
                'x': np.arange(self.start, self.end),
                'y': total_coverage,
                'height': total_coverage
            })
            
            p.vbar(x='x', top='y', width=1, source=source,
                fill_color='gray', line_color=None)
            
            p.xaxis.visible = False
            p.yaxis.visible = True
            min_val = 0
            max_val = int(max(total_coverage))
            p.yaxis.ticker = [min_val, max_val]
    
        p.grid.visible = False

        p.toolbar_location = None
            
        # self.plots.append(p)
        # self._update_layout()
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
        p.tags = [{"type": "alignment",
                   "bam_path": bam_path,
                    }]

        if (self.end - self.start) >= 10000:
            p.text(x=[self.start + (self.end - self.start)/2],
                    y=[25],
                    text=['Zoom in to view alignments'],
                    text_font_size='8pt',
                    text_baseline="middle",
                    text_align="center",
                    text_color='gray')
        else:
            
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
        p.tags = [{"type": "annotations"}]
        
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
    
    def serve(self, port=5006, num_procs=1, show=True):
        """
        Serve the Pileupy browser as a Bokeh application.
        
        :param port: Port to serve the application on.
        :param num_procs: Number of worker processes.
        :param show: Whether to open the browser automatically.
        """
        pass
        def bkapp(doc):
            # Create a fresh layout for each document
            doc_layout = column(self.plots)
            doc.add_root(doc_layout)

        server = Server({'/': bkapp}, port=port, num_procs=num_procs)
        server.start()

        if show:
            print(f"Opening Pileupy browser at http://localhost:{port}/")
            server.io_loop.add_callback(server.show, "/")

        server.io_loop.start()

def _chr_sort_key(chrom):
    """Custom sorting key for chromosome names"""
    # Remove 'chr' prefix if it exists
    chrom = chrom.replace('chr', '')
    
    # Handle numbered chromosomes
    if chrom.isdigit():
        return (0, int(chrom))  # (priority, number)
    # Handle X chromosome
    elif chrom == 'X':
        return (1, 23)  # After numbered chromosomes
    # Handle Y chromosome
    elif chrom == 'Y':
        return (1, 24)  # After X
    # Handle M/MT (mitochondrial)
    elif chrom in ['M', 'MT']:
        return (1, 25)
    # Handle other contigs/scaffolds
    else:
        return (2, chrom)  # Put other contigs last, sorted alphabetically