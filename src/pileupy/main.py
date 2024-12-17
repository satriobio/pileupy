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


class TrackControl:
    def __init__(self, chr_lengths, chrom, start, end, update_callback):
        self.update_callback = update_callback        
        self._fetch(chr_lengths, chrom, start, end)
        self._figure()
    
    def _fetch(self, chr_lengths, chrom, start, end):
        self.chr_lengths = chr_lengths        
        self.chroms = sorted(list(chr_lengths.keys()), key=_chr_sort_key)
        self.chrom = chrom
        self.start = start
        self.end = end
    
    def _update(self):
        # This function can be used internally
        self.update_callback(self.chrom, self.start, self.end)  # Call the callback when needed
    
    def _figure(self):
        chrom_select = Select(
            value=self.chrom,
            options=self.chroms,
            width=100
        )

        start_input = TextInput(
            value=str(self.start),
            width=120
        )
        
        end_input = TextInput(
            value=str(self.end),
            width=120
        )

        zoom_out = Button(label="-", button_type="primary", width=15)
        zoom_in = Button(label="+", button_type="primary", width=15)

        # Define callbacks
        def update_region(attr, old, new):
            """Update the viewed region"""
            new_start = int(start_input.value.replace(',', ''))
            new_end = int(end_input.value.replace(',', ''))                    
            
            self.start = new_start
            self.end = new_end

            self._update()

        def update_chrom(attr, old, new):
            """Handle chromosome selection"""
            self.chrom = new
            self.start = 0
            self.end = self.chr_lengths.get(new, 0)
            start_input.value = str(self.start)
            end_input.value = str(self.end)

            self._update()

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
            
            self.start = new_start
            self.end = new_end
            
            start_input.value = str(new_start)
            end_input.value = str(new_end)

            self._update()

        # Connect callbacks
        chrom_select.on_change('value', update_chrom)
        start_input.on_change('value', update_region)
        end_input.on_change('value', update_region)
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
            margin=(0, 0, 0, 20),
            sizing_mode="stretch_width",  # Make controls stretch horizontally
            width_policy="max"
        )

        self.figure = column(
            controls, 
            sizing_mode="stretch_width",
            width_policy="max",
            # height_policy="fit"
        )

    
    def _renderers(self):
        pass

class TrackIdiogram:
    def __init__(self, cytobands, chr_lengths, chrom, start, end, shared_x_range):
        self.stain_colors = {
            'gneg': '#ffffff',
            'gpos25': '#c0c0c0',
            'gpos50': '#808080',
            'gpos75': '#404040',
            'gpos100': '#000000',
            'acen': '#963232',
            'gvar': '#000000',
            'stalk': '#963232',
        }
        
        # Initialize data sources
        self.source_backbone = ColumnDataSource({
            'x': [], 'y': [], 'width': [], 'height': []
        })
        self.source_bands = ColumnDataSource({
            'x': [], 'y': [], 'width': [], 'height': [], 'color': []
        })
        self.source_centromeres = ColumnDataSource({
            'xs': [], 'ys': [], 'color': []
        })
        self.source_view = ColumnDataSource({
            'x': [], 'y': [], 'width': [], 'height': []
        })
        
        self._figure()
        self._renderers()
        self._fetch(cytobands, chr_lengths, chrom, start, end)

    def _fetch(self, cytobands, chr_lengths, chrom, start, end):
        """Prepare data for plotting"""
        self.cytobands = cytobands
        self.chr_lengths = chr_lengths
        self.chrom = chrom

        self.start = start
        self.end = end
        
        # Filter cytobands for current chromosome
        self.chr_cytobands = [c for c in self.cytobands if c['chrom'] == self.chrom]
        self.chr_length = self.chr_lengths.get(self.chrom, 0)
        
        # Prepare backbone data
        backbone_data = {
            'x': [self.chr_length/2],
            'y': [15],
            'width': [self.chr_length],
            'height': [10]
        }
        self.source_backbone.data = backbone_data
        
        # Prepare band data
        band_data = {
            'x': [], 'y': [], 'width': [], 'height': [], 'color': []
        }
        centromere_data = {
            'xs': [], 'ys': [], 'color': []
        }
        
        for band in self.chr_cytobands:
            width = band['end'] - band['start']
            x = band['start'] + width/2
            
            if band['stain'] == 'acen':
                # Centromere data
                if band['name'].startswith('p'):
                    xs = [band['start'], band['end'], band['start']]
                    ys = [10, 15, 20]
                else:
                    xs = [band['start'], band['end'], band['end']]
                    ys = [15, 20, 10]
                centromere_data['xs'].append(xs)
                centromere_data['ys'].append(ys)
                centromere_data['color'].append(self.stain_colors['acen'])
            else:
                # Regular band data
                band_data['x'].append(x)
                band_data['y'].append(15)
                band_data['width'].append(width)
                band_data['height'].append(10)
                band_data['color'].append(self.stain_colors.get(band['stain'], '#ffffff'))
        
        self.source_bands.data = band_data
        self.source_centromeres.data = centromere_data
        
        # Prepare view indicator data
        view_data = {
            'x': [(self.start + self.end) / 2],
            'y': [15],
            'width': [self.end - self.start],
            'height': [10]
        }
        self.source_view.data = view_data

    def _update(self, ref, chrom, start, end, shared_x_range):
        """Update the visualization with new data"""
        
        # Update the shared x_range
        self._fetch(self.cytobands, self.chr_lengths, chrom, start, end)

        self.figure.x_range.start = 0
        self.figure.x_range.end = self.chr_length

    def _figure(self):
        """Create and configure the figure"""
        fig = bplt.figure(
            width_policy="max",  # Allow figure to stretch horizontally
            height=30,  # Fixed height
            sizing_mode="stretch_width",  # Make figure responsive to container width
            tools=''
        )
        
        # Configure axes
        fig.xaxis.visible = False
        fig.yaxis.visible = False
        fig.grid.visible = False
        
        # Disable all interactions
        fig.toolbar_location = None
        fig.tools = []
        fig.toolbar.tools = []
        
        self.figure = fig
        return fig

    def _renderers(self):
        """Set up all renderers using data sources"""
        
        # Draw chromosome backbone
        self.figure.rect(
            x='x', y='y',
            width='width', height='height',
            fill_color='white',
            line_color='black',
            source=self.source_backbone
        )
        
        # Draw regular bands
        self.figure.rect(
            x='x', y='y',
            width='width', height='height',
            fill_color='color',
            line_color='black',
            line_width=0.5,
            source=self.source_bands
        )
        
        # Draw centromeres
        self.figure.patches(
            xs='xs', ys='ys',
            fill_color='color',
            line_color='black',
            source=self.source_centromeres
        )
        
        # Draw current view indicator
        self.figure.rect(
            x='x', y='y',
            width='width', height='height',
            fill_color='red',
            fill_alpha=0.2,
            line_color='red',
            line_width=1,
            source=self.source_view
        )

class TrackReference:
    def __init__(self, ref, chrom, start, end, shared_x_range):
        self.base_colors = {
            'A': '#00FF00',  # Green
            'T': '#FF0000',  # Red
            'C': '#0000FF',  # Blue
            'G': '#FFA500'   # Orange
        }
        self.source_block = ColumnDataSource({
            'x': [], 
            'y': [],
            'text': [],
            'colors': []
        })
        self.source_letter = ColumnDataSource({
            'x': [], 
            'y': [],
            'text': [],
            'colors': []
        })
        self.source_message = ColumnDataSource({
            'x': [],
            'y': [],
            'text': []
        })
        
        self.shared_x_range = shared_x_range
        self._figure(shared_x_range)
        self._renderers()
        self._fetch(ref, chrom, start, end)

    def _fetch(self, ref, chrom, start, end):
        """Prepare data for plotting by updating the sources"""
        self.start = start
        self.end = end
        span = self.end - self.start

        x = np.arange(start, end)
        y = [2] * len(x)
        
        # Clear all sources
        empty_data = {
            'x': [], 
            'y': [],
            'text': [],
            'colors': []
        }
        self.source_block.data = empty_data.copy()
        self.source_letter.data = empty_data.copy()
        self.source_message.data = {'x': [], 'y': [], 'text': []}

        if span >= 10000:
            # Show "zoom in" message
            self.source_message.data = {
                'x': [start + span/2],
                'y': [2],
                'text': ['Zoom in to see reference sequence']
            }
        else:
            # Fetch sequence data
            self.ref_seq = ref.fetch(chrom, start, end)
            sequence = list(self.ref_seq.upper())
            colors = [self.base_colors.get(base.upper(), 'gray') for base in self.ref_seq]

            if span >= 70:  # Show blocks
                self.source_block.data = {
                    'x': x,
                    'y': y,
                    'text': sequence,
                    'colors': colors
                }
            else:  # Show letters
                self.source_letter.data = {
                    'x': x,
                    'y': y,
                    'text': sequence,
                    'colors': colors
                }

    def _update(self, ref, chrom, start, end, shared_x_range):
        """Update the visualization with new data"""

        # print(f"Updating track reference with {chrom}:{start}-{end}, {shared_x_range}")
        self._fetch(ref, chrom, start, end)
        
        # Update the shared x_range
        self.shared_x_range.start = start
        self.shared_x_range.end = end

    def _figure(self, shared_x_range):
        """Create and configure the figure"""
        fig = bplt.figure(width_policy="max", height=50, tools='', x_range=shared_x_range)
        
        # Configure axes
        fig.xaxis.visible = True
        fig.xaxis.axis_label_text_font_size = '10pt'
        fig.xaxis[0].formatter.use_scientific = False
        fig.xaxis.axis_line_color = 'black'
        fig.xaxis.major_tick_line_color = 'black'
        fig.xaxis.minor_tick_line_color = None
        fig.above = fig.below  # Move axis to top
        fig.below = []  # Remove bottom axis

        fig.yaxis.visible = False
        fig.grid.visible = False
    
        # Disable all interactions
        fig.toolbar_location = None
        fig.tools = []
        fig.toolbar.tools = []

        self.figure = fig
        return fig

    def _renderers(self):
        """Set up all renderers once during initialization"""

        # Add message renderer
        self.figure.text(
            x='x',
            y='y',
            text='text',
            text_font_size='8pt',
            text_baseline="middle",
            text_align="center",
            text_color='gray',
            source=self.source_message
        )
        
        # Add block renderer
        self.figure.rect(
            x='x', 
            y='y', 
            width=1, 
            height=2,
            fill_color='colors', 
            line_color=None, 
            source=self.source_block
        )
        
        # Add letter renderer
        self.figure.text(
            x='x', 
            y='y', 
            text='text', 
            text_color='colors', 
            source=self.source_letter,
            text_baseline="middle", 
            text_align="center"
        )

class TrackAlignment:
    def __init__(self, ref, bam, chrom, start, end, shared_x_range):
        # Initialize data sources
        self.base_colors = {
            'A': '#00FF00',  # Green
            'T': '#FF0000',  # Red
            'C': '#0000FF',  # Blue
            'G': '#FFA500'   # Orange
        }

        self.source_coverage = ColumnDataSource({
            'x': [],
            'y': [],
            'height': []
        })

        self.source_reads = ColumnDataSource({
            'xs': [], 'ys': [], 'colors': [],
            'read_names': [], 'mapping_qualities': [], 
            'cigar_strings': [], 'read_lengths': []
        })
        
        self.source_snps_rect = ColumnDataSource({
            'x': [], 'y': [], 'height': [], 'color': []
        })
        
        self.source_snps_text = ColumnDataSource({
            'x': [], 'y': [], 'color': [], 'text': []
        })
        
        self.source_message = ColumnDataSource({
            'x': [], 'y': [], 'text': []
        })
        
        # Set fixed dimensions for reads
        self.read_height = 10  # Fixed read height
        self.spacing = 6      # Fixed spacing between reads
        self.default_view_height = 400  # Default visible height
        
        self.shared_x_range = shared_x_range
        self._figure(shared_x_range)
        self._renderers()
        self._fetch(ref, bam, chrom, start, end)

    def _fetch(self, ref, bam, chrom, start, end):
        self.bam = bam
        self.start = start
        self.end = end

        self.ref_seq = ref.fetch(chrom, start, end)

        if (self.end - self.start) >= 10000:
            self.source_coverage.data = {
                'x': [],
                'y': [],
                'height': []
            }
        else:
            coverage = bam.count_coverage(chrom, start, end)
            total_coverage = np.sum(coverage, axis=0)  # Sum across bases
            self.source_coverage.data = {
                'x': np.arange(self.start, self.end),
                'y': total_coverage,
                'height': total_coverage
            }

        # Prepare data containers
        message_data = {
            'x': [], 'y': [], 'text': []
        }
        read_data = {
            'xs': [], 'ys': [], 'colors': [],
            'read_names': [], 'mapping_qualities': [],
            'cigar_strings': [], 'read_lengths': []
        }
        snp_rect_data = {
            'x': [], 'y': [], 'height': [], 'color': []
        }
        snp_text_data = {
            'x': [], 'y': [], 'color': [], 'text': []
        }

        if (self.end - self.start) >= 10000:
            self.source_message.data = {
                'x': [start + (end - start)/2],
                'y': [self.default_view_height/2],
                'text': ['Zoom in to view alignments']
            }
            self.source_reads.data = read_data
            self.source_snps_rect.data = snp_rect_data
            self.source_snps_text.data = snp_text_data
            return
        else:
            self.reads = list(bam.fetch(chrom, start, end))
            self.source_message.data = message_data

        # Track occupied positions for read placement
        occupied_positions = []  # List of (start, end) tuples for each row

        # Process reads and collect data
        for read in self.reads:
            if not read.is_unmapped and read.reference_length:
                read_start = read.reference_start
                read_end = read.reference_start + read.reference_length
                
                # Find first available row for read placement
                row = 0
                placed = False
                
                # Try to place read in existing rows
                for i, positions in enumerate(occupied_positions):
                    can_place = True
                    for pos_start, pos_end in positions:
                        if not (read_end < pos_start or read_start > pos_end):
                            can_place = False
                            break
                    
                    if can_place:
                        row = i
                        occupied_positions[i].append((read_start, read_end))
                        placed = True
                        break
                
                # If no existing row works, create new row
                if not placed:
                    row = len(occupied_positions)
                    occupied_positions.append([(read_start, read_end)])

                # Calculate y position (from top down)
                y_pos = self.default_view_height - (row * (self.read_height + self.spacing))
                
                # Calculate arrow points
                x_start = read.reference_start - 1
                x_end = read.reference_start + read.reference_length - 1
                y_top = y_pos + self.read_height
                y_bottom = y_pos

                # Make arrow_width proportional to read height
                if abs(start-end) > 4000:    
                    arrow_width = 0
                else:
                    arrow_width = abs(start-end)/150  # This makes the arrow width equal to the height for a 45-degree angle
                y_mid = y_pos + (self.read_height / 2)  # Calculate middle y-position

                if read.is_reverse:
                    # Left-pointing arrow (reverse strand)
                    xs = [
                        x_start, x_end,          # Bottom line
                        x_end, x_end,            # Right edge
                        x_end, x_start + arrow_width,  # Top line
                        x_start + arrow_width/2, x_start,  # Arrow point top
                        x_start, x_start + arrow_width/2,  # Arrow point bottom
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
                        x_end - arrow_width/2, x_end,    # Arrow point bottom
                        x_end, x_end - arrow_width/2,    # Arrow point top
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
                read_data['colors'].append('#ABABAB')
                
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
                                        if abs(self.end - self.start) < 70:
                                            # Text SNPs (when zoomed in)
                                            snp_text_data['x'].append(read.reference_start + i)
                                            snp_text_data['y'].append(y_pos + self.read_height/2)  # Center vertically in read
                                            snp_text_data['color'].append(self.base_colors.get(query_base, '#808080'))
                                            snp_text_data['text'].append(query_base)
                                        else:
                                            # Rectangle SNPs (when zoomed out)
                                            snp_rect_data['x'].append(read.reference_start + i)
                                            snp_rect_data['y'].append(y_top - self.read_height/2)  # Use exact same y_pos as read
                                            snp_rect_data['height'].append(self.read_height)  # Use exact same height as read
                                            snp_rect_data['color'].append(self.base_colors.get(query_base, '#808080'))
                                    
                            ref_pos += length  # Update reference position after the block
                            query_pos += length
                        elif op == 1:  # Insertion
                            query_pos += length
                        elif op == 2:  # Deletion
                            ref_pos += length  # Update reference position for deletions
                            
        for read in self.reads:
            if not read.is_unmapped and read.reference_length:
                read_data['read_names'].append(read.query_name)
                read_data['mapping_qualities'].append(read.mapping_quality)
                read_data['cigar_strings'].append(read.cigarstring)
                read_data['read_lengths'].append(len(read.query_sequence) if read.query_sequence else 0)
        
        # Update data sources
        self.source_reads.data = read_data
        self.source_snps_rect.data = snp_rect_data
        self.source_snps_text.data = snp_text_data

        min_y = min([min(y) for y in read_data['ys']]) - self.read_height
        self.figure_aln.y_range = Range1d(0, self.default_view_height, bounds=(min_y, self.default_view_height))

    def _update(self, ref, chrom, start, end, shared_x_range):
        """Update the visualization with new data"""
        self._fetch(ref, self.bam, chrom, start, end)
        
        # Update the shared x_range
        self.shared_x_range.start = start
        self.shared_x_range.end = end

    def _figure(self, shared_x_range):
        
        fig_cov = bplt.figure(
            height=50,
            tools="",
            x_range=shared_x_range,
            sizing_mode="stretch_width",
            width_policy="max"
        )

        fig_cov.xaxis.visible = False
        fig_cov.yaxis.visible = False
        fig_cov.grid.visible = False
        fig_cov.toolbar.logo = None

        fig_aln = bplt.figure(
            height=self.default_view_height,
            tools="ywheel_pan",
            active_scroll="ywheel_pan",
            x_range=shared_x_range,
            sizing_mode="stretch_width",
            width_policy="max"
        )
        
        fig_aln.xaxis.visible = False
        fig_aln.yaxis.visible = False
        fig_aln.grid.visible = False
        fig_aln.toolbar.logo = None
        

        self.figure_aln = fig_aln
        self.figure_cov = fig_cov
        self.figure = column(
            self.figure_cov,
            self.figure_aln, 
            sizing_mode="stretch_width",
            width_policy="max"
        )
        return self.figure

    def _renderers(self):
        self.figure_cov.vbar(
            x='x',
            top='y',
            width=1,
            source=self.source_coverage,
            fill_color='gray',
            line_color=None
        )

        # Add message renderer
        self.figure_aln.text(
            x='x', y='y', text='text',
            text_font_size='8pt',
            text_baseline="middle",
            text_align="center",
            text_color='gray',
            source=self.source_message
        )
        
        # Add read patches renderer
        patches = self.figure_aln.patches(
            xs='xs', ys='ys',
            fill_color='colors',
            line_color=None,
            alpha=0.8,
            source=self.source_reads
        )
        
        # Add SNP renderers (both rect and text)
        self.figure_aln.rect(
            x='x',
            y='y',
            width=1,
            height='height',
            color='color',
            line_color=None,
            alpha=0.8,
            source=self.source_snps_rect
        )
        
        self.figure_aln.text(
            x='x',
            y='y',
            text='text',
            text_font_size='10pt',
            text_font_style='bold',
            text_baseline="middle",
            text_align="center",
            text_color='color',
            source=self.source_snps_text
        )
        
        # Add hover tool
        hover = HoverTool(
            tooltips=[
                ('Read Name', '@read_names'),
                ('Mapping Quality', '@mapping_qualities'),
                ('CIGAR', '@cigar_strings'),
                ('Length', '@read_lengths')
            ],
            renderers=[patches]
        )
        self.figure_aln.add_tools(hover)
        self.figure_aln.add_tools(TapTool())

class TrackAnnotation:
    def __init__(self):
        pass

class TrackContinuous:
    def __init__(self):
        pass

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
        self.tracks = []
        # Create initial tracks based on input type
        if genome:
            self._create_base_tracks()
        else:
            self.add_track_reference()
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

        # self.track_control()
        # self.track_idiogram()
        self.add_track_control()
        self.add_track_idiogram()
        self.add_track_reference()
        # self.track_genes()
        # self._update_layout()
        # self._update_tracks()

    def add_track_control(self):
        track_control = TrackControl(self.chr_lengths, self.chrom, self.start, self.end, self._update_tracks)
        # self.chrom = track_control.chrom
        # self.start = track_control.start
        # self.end = track_control.end
        self.tracks.append(track_control)

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
            # addons = self.plots[4:]

            
            # Clear existing plots and recreate base tracks
            # self.plots = []
            self._create_base_tracks()
            
            # Process addon tracks in batch
            # if addons:
            #     for addon in addons:
            #         if len(addon.tags) > 0 and addon.tags[0].get('type') == 'alignment':
            #             self.track_alignment(addon.tags[0].get('bam_path'))
            #         if len(addon.tags) > 0 and addon.tags[0].get('type') == 'annotation':
            #             self.track_alignment(addon.tags[0].get('bed_path'))
            # self._update_layout()

            # if curdoc().session_context:
            #     curdoc().clear()
            #     curdoc().add_root(self.layout)

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
            margin=(0, 0, 0, 20),
            sizing_mode="stretch_width",  # Make controls stretch horizontally
            width_policy="max"
        )

        # self.plots.append(controls_layout)        
        # self._update_layout()
        controls_layout = column(controls, sizing_mode="stretch_width")
        self.tracks.append(controls_layout) 
        # return controls_layout
    
    def add_track_idiogram(self):
        """Create chromosome idiogram track with cytobands"""
        track_idiogram = TrackIdiogram(self.cytobands, self.chr_lengths, self.chrom, self.start, self.end, self.shared_x_range)
        self.tracks.append(track_idiogram)
        # Create the figure first with shared x range
        # p = bplt.figure(width=800, height=50, tools='', x_range=self.shared_x_range)
        
        # # Filter cytobands for current chromosome
        # chr_cytobands = [c for c in self.cytobands if c['chrom'] == self.chrom]
        # chr_length = self.chr_lengths.get(self.chrom, 0)
        
        # if not chr_cytobands or chr_length == 0:
        #     print("No cytoband data found for chromosome", self.chrom)
        #     return p
        
        # stain_colors = {
        #     'gneg': '#ffffff',
        #     'gpos25': '#c0c0c0',
        #     'gpos50': '#808080',
        #     'gpos75': '#404040',
        #     'gpos100': '#000000',
        #     'acen': '#963232',
        #     'gvar': '#000000',
        #     'stalk': '#963232',
        # }
        
        # # Set up the plot
        # p.height = 50
        # p.y_range.start = 0
        # p.y_range.end = 30
        
        # # Fix x_range to show full chromosome, regardless of current view
        # p.x_range = Range1d(0, chr_length)
        
        # # Draw chromosome backbone
        # p.rect(x=chr_length/2,
        #        y=15,
        #        width=chr_length,
        #        height=10,
        #        fill_color='white',
        #        line_color='black')
        
        # # Draw cytobands
        # for band in chr_cytobands:
        #     width = band['end'] - band['start']
        #     x = band['start'] + width/2
            
        #     if band['stain'] == 'acen':
        #         # Draw centromere as triangle
        #         if band['name'].startswith('p'):
        #             xs = [band['start'], band['end'], band['start']]
        #             ys = [10, 15, 20]
        #         else:
        #             xs = [band['start'], band['end'], band['end']]
        #             ys = [15, 20, 10]
                
        #         p.patches(
        #             xs=[xs],
        #             ys=[ys],
        #             fill_color=stain_colors['acen'],
        #             line_color='black'
        #         )
        #     else:
        #         # Draw regular band
        #         p.rect(
        #             x=x,
        #             y=15,
        #             width=width,
        #             height=10,
        #             fill_color=stain_colors.get(band['stain'], '#ffffff'),
        #             line_color='black',
        #             line_width=0.5
        #         )
        
        # # Add current view indicator (red line)
        # p.rect(
        #     x=(self.start + self.end) / 2,  # center x position
        #     y=15,  # center y position (same as chromosome backbone)
        #     width=self.end - self.start,  # width spans the current view
        #     height=10,  # same height as chromosome backbone
        #     fill_color='red',
        #     fill_alpha=0.2,
        #     line_color='red',
        #     line_width=1
        # )
        
        # # Configure axes
        # p.xaxis.visible = False
        # p.yaxis.visible = False
        # p.grid.visible = False
        
        # # Disable all interactions
        # p.toolbar_location = None
        # p.tools = []
        # p.toolbar.tools = []
        
        # self.plots.append(p)
        # self._update_layout()
        # return p
    
    def add_track_reference(self):
        """Create the reference sequence track"""
        # Create the figure first with shared x range
        track_reference = TrackReference(self.ref, self.chrom, self.start, self.end, self.shared_x_range)
        self.tracks.append(track_reference)
        # self._update_layout()
        # return track_reference
        
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
        # track_reference = TrackReference(self.ref, self.chrom, self.start, self.end, self.shared_x_range)
        # self.tracks.append(track_reference)
        bam = pysam.AlignmentFile(bam_path)
        track_alignment = TrackAlignment(self.ref, bam, self.chrom, self.start, self.end, self.shared_x_range)
        self.tracks.append(track_alignment)

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
            self.layout = column(
                self.plots,
                sizing_mode="stretch_width",  # Make layout stretch to container width
                width_policy="max",  # Allow layout to take maximum width
                height_policy="fit"  # Fit height to content
            )
            print(f"Layout updated with {len(self.plots)} plots")
        except Exception as e:
            print(f"Error updating layout: {str(e)}")

    def _update_tracks(self, chrom, start, end):
        """Update the tracks after adding new tracks"""
        # Update the shared x_range
        self.shared_x_range.start = start
        self.shared_x_range.end = end
        
        for track in self.tracks[1:]:
            track._update(self.ref, chrom, start, end, self.shared_x_range)
        # layout = column([track.plot for track in self.tracks])
        # self.plots = []
        # self._create_base_tracks()
        # self._update_layout()
        # pass

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
        def bkapp(doc):
            # Create a fresh layout for each document
            doc_layout = column(
                [track.figure for track in self.tracks],
                sizing_mode="stretch_both",  # Make layout responsive in both directions
                width_policy="max",
                height_policy="max"
            )
            doc.add_root(doc_layout)

            # Add responsive sizing to document
            # doc.theme = bplt.document.Theme(
            #     json={
            #         'attrs': {
            #             'Figure': {
            #                 'sizing_mode': 'stretch_width',
            #                 'width_policy': 'max',
            #                 'height_policy': 'fit'
            #             },
            #             'Column': {
            #                 'sizing_mode': 'stretch_width',
            #                 'width_policy': 'max',
            #                 'height_policy': 'fit'
            #             },
            #             'Row': {
            #                 'sizing_mode': 'stretch_width',
            #                 'width_policy': 'max',
            #                 'height_policy': 'fit'
            #             }
            #         }
            #     }
            # )

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