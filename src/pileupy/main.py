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
        self.filter = []        
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
        self.update_callback(self.chrom, self.start, self.end, self.filter)  # Call the callback when needed
    
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

class TrackGene:
    def __init__(self, genes, chrom, start, end, shared_x_range):
        # Store initial parameters
        self.genes = genes
        self.chrom = chrom
        self.start = start
        self.end = end
        
        # Initialize data sources
        self.source_lines = ColumnDataSource({
            'x': [], 'y': []
        })
        
        self.source_exons = ColumnDataSource({
            'x': [], 'y': [], 'width': [], 'height': []
        })

        self.source_labels = ColumnDataSource({
            'x': [], 'y': [], 'text': []
        })

        self.shared_x_range = shared_x_range
        self._figure(shared_x_range)
        self._renderers()
        self._fetch(genes, chrom, start, end)
    
    def _fetch(self, genes, chrom, start, end):
        """Prepare data for plotting"""
        self.start = start
        self.end = end
        self.genes = genes
        
        visible_genes = [g for g in genes 
                        if g['chrom'] == chrom and 
                        not (g['end'] < start or g['start'] > end)]

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

        for gene in visible_genes:
            # Add gene body line coordinates if zoomed in enough
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
        
        # Update data sources
        # self.source_lines.data = {
        #     'x': gene_lines_x,
        #     'y': gene_lines_y
        # }

        self.source_exons.data = {
            'x': exon_xs,
            'y': exon_ys,
            'width': exon_widths,
            'height': [exon_height] * len(exon_xs)
        }

        self.source_labels.data = {
            'x': text_xs,
            'y': text_ys,
            'text': text_labels
        }

    def _update(self, ref, chrom, start, end, shared_x_range):
        """Update the visualization with new data"""
        self._fetch(self.genes, chrom, start, end)
        
        # Update the shared x_range
        self.shared_x_range.start = start
        self.shared_x_range.end = end

    def _figure(self, shared_x_range):
        fig = bplt.figure(width_policy="max", height=50, tools='', x_range=shared_x_range)
        fig.y_range = Range1d(0, 50)

        fig.xaxis.visible = False
        fig.yaxis.visible = False
        fig.grid.visible = False

        fig.toolbar_location = None
        fig.tools = []
        fig.toolbar.tools = []

        self.figure = fig
        return fig

    def _renderers(self):
        # Add gene body lines
        self.figure.line(
            x='x', 
            y='y', 
            source=self.source_lines, 
            line_color='blue', 
            line_width=1
        )

        # Add exons
        self.figure.rect(
            x='x', 
            y='y', 
            width='width', 
            height='height',
            color='blue', 
            line_color=None, 
            source=self.source_exons
        )

        # Add gene labels
        self.figure.text(
            x='x', 
            y='y', 
            text='text',
            text_font_size='8pt', 
            text_baseline="middle", 
            text_align="center",
            source=self.source_labels
        )

class TrackAlignment:
    def __init__(self, ref, bam, chrom, start, end, shared_x_range, height=400, update_callback=None):
        # Initialize data sources
        self.update_callback = update_callback
        self.filter = []  
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
        self.default_view_height = height  # Default visible height
        
        self.shared_x_range = shared_x_range
        self._figure(shared_x_range)
        self._renderers()
        self._fetch(ref, bam, chrom, start, end)

    def _fetch(self, ref, bam, chrom, start, end):
        self.bam = bam
        self.chrom = chrom
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
                        if op == 0 or op == 7 or op == 8:  # M, = or X (match/mismatch)
                            for i in range(length):
                                if (ref_pos + i >= 0 and 
                                    ref_pos + i < len(self.ref_seq) and 
                                    query_pos + i < len(read.query_sequence)):

                                    ref_base = self.ref_seq[ref_pos + i].upper()
                                    query_base = read.query_sequence[query_pos + i].upper()
                                        
                                    if ref_base != query_base:
                                        if abs(self.end - self.start) < 70:
                                            # Text SNPs (when zoomed in)
                                            snp_text_data['x'].append(read.reference_start + query_pos + i)
                                            snp_text_data['y'].append(y_pos + self.read_height/2)
                                            snp_text_data['color'].append(self.base_colors.get(query_base, '#808080'))
                                            snp_text_data['text'].append(query_base)
                                        else:
                                            # Rectangle SNPs (when zoomed out)
                                            snp_rect_data['x'].append(read.reference_start + query_pos + i)
                                            snp_rect_data['y'].append(y_top - self.read_height/2)
                                            snp_rect_data['height'].append(self.read_height)
                                            snp_rect_data['color'].append(self.base_colors.get(query_base, '#808080'))
                            ref_pos += length
                            query_pos += length
                        elif op == 1:  # I (Insertion)
                            snp_text_data['x'].append(read.reference_start + ref_pos)
                            snp_text_data['y'].append(y_pos + self.read_height/2)
                            snp_text_data['color'].append('#800080')  # Purple
                            snp_text_data['text'].append('I')
                            query_pos += length
                        elif op == 2:  # D (Deletion)
                            # Add light gray rectangle for deletion
                            snp_rect_data['x'].append(read.reference_start + ref_pos)
                            snp_rect_data['y'].append(y_top - self.read_height/2)
                            snp_rect_data['height'].append(self.read_height)
                            snp_rect_data['color'].append('#D3D3D3')
                            ref_pos += length
                        elif op == 4:  # S (Soft clipping)
                            query_pos += length
                        elif op == 5:  # H (Hard clipping)
                            pass
                        elif op == 3:  # N (Skipped region from reference)
                            ref_pos += length
            # break
        
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
        self.figure_aln.y_range = Range1d(0, self.default_view_height + 20, bounds=(min_y, self.default_view_height + 20))

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

        self.source_reads.selected.on_change('indices', lambda attr, old, new: self.on_read_click(attr, old, new, self.source_reads))
        

        self.figure_aln.add_tools(hover)
        self.figure_aln.add_tools(TapTool())
    
    def on_read_click(self, attr, old, new, source):
        if new:  # If there are selected indices
            reads = []
            for x in new:
                reads.append(source.data['read_names'][x])
                print(f"new {source.data['read_names'][x]}")
            
            self.update_callback(self.chrom, self.start, self.end, reads)

class TrackAnnotation:
    def __init__(self, bed_path, chrom, start, end, shared_x_range):
        # Store initial parameters
        self.bed_path = bed_path
        self.chrom = chrom
        self.start = start
        self.end = end
        
        # Initialize data source
        self.source_rect = ColumnDataSource({
            'x': [], 'y': [], 'width': [], 'height': [], 'colors': []
        })

        self.source_text = ColumnDataSource({
            'x': [], 'y': [], 'text': [], 'colors': []
        })

        self.shared_x_range = shared_x_range
        self._figure(shared_x_range)
        self._renderers()
        self._fetch(bed_path, chrom, start, end)

    def _fetch(self, bed_path, chrom, start, end):
        """Prepare data for plotting"""
        self.start = start
        self.end = end
        
        # Load and filter annotations
        annotations = load_bed(bed_path)
        visible_annotations = [f for f in annotations 
                             if f['chrom'] == chrom and 
                             not (f['end'] < start or f['start'] > end)]
        
        if not visible_annotations:
            # Show "no data" message
            # self.source_text.data = {
            #     'x': [start + (end - start)/2],
            #     'y': [15],
            #     'text': ['No annotations in this region'],
            #     'colors': ['gray']
            # }
            self.source_rect.data = {
                'x': [], 'y': [], 'width': [], 'height': [], 'colors': []
            }
        else:
            # Prepare data for rectangles
            self.source_rect.data = {
                'x': [(ann['start'] + ann['end'])/2 for ann in visible_annotations],
                'y': [35] * len(visible_annotations),  # Fixed y position at 35
                'width': [ann['end'] - ann['start'] for ann in visible_annotations],
                'height': [20] * len(visible_annotations),  # Fixed height of 20
                'colors': ['red'] * len(visible_annotations)
            }
            
            # Prepare data for text labels
            self.source_text.data = {
                'x': [(ann['start'] + ann['end'])/2 for ann in visible_annotations],
                'y': [15] * len(visible_annotations),  # Fixed y position at 15
                'text': [ann['name'] for ann in visible_annotations],
                'colors': ['black'] * len(visible_annotations)
            }

    def _update(self, ref, chrom, start, end, shared_x_range):
        """Update the visualization with new data"""
        self._fetch(self.bed_path, chrom, start, end)
        
        # Update the shared x_range
        self.shared_x_range.start = start
        self.shared_x_range.end = end
    
    def _figure(self, shared_x_range):
        """Create and configure the figure"""
        fig = bplt.figure(width_policy="max", height=50, tools='', x_range=shared_x_range)
        fig.y_range = Range1d(0, 50)

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
        # Add rectangle renderer for annotation blocks
        self.figure.rect(
            x='x',
            y='y',
            width='width',
            height='height',
            color='colors',
            alpha=0.7,
            source=self.source_rect
        )
        
        # Add text renderer for annotation names
        self.figure.text(
            x='x',
            y='y',
            text='text',
            text_font_size='8pt',
            text_baseline="middle",
            text_align="center",
            source=self.source_text
        )

class TrackDataFrame:
    def __init__(self, df, chrom, start, end, shared_x_range, height=100):
        self.height = height
        # Store initial parameters
        self.df = df  # Store df as instance variable
        self.chrom = chrom
        self.start = start
        self.end = end
        
        # Initialize data source
        self.source_lines = ColumnDataSource({
            'xs': [], 'ys': [], 'colors': [], 'read_names': []
        })

        self.shared_x_range = shared_x_range
        self._figure(shared_x_range)
        self._renderers()
        self._fetch(df, chrom, start, end)

    def _fetch(self, df, chrom, start, end):
        """Prepare data for plotting"""
        self.start = start
        self.end = end
        
        # Filter data for current chromosome and view range
        mask = (df['chrom'] == chrom) & \
               (df['start'] <= end) & \
               (df['end'] >= start)
        visible_data = df[mask].copy()
        
        if visible_data.empty:
            self.source_lines.data = {
                'xs': [], 'ys': [], 'colors': [], 'read_names': []
            }
        else:
            # Group data by 'read_names' for multi-line plotting
            xs = []
            ys = []
            colors = []
            read_names = []

            for read_name, group in visible_data.groupby('read_names'):
                xs.append(group['start'].tolist())  # x-coordinates for this line
                ys.append(group['value'].tolist())  # y-coordinates for this line
                colors.append('blue')  # Default color
                read_names.append(read_name)  # Store the read name for each line

            self.source_lines.data = {
                'xs': xs,
                'ys': ys,
                'colors': colors,
                'read_names': read_names
            }

            # print(self.source_lines.data)

    def _update(self, ref, chrom, start, end, shared_x_range):
        """Update the visualization with new data"""
        self._fetch(self.df, chrom, start, end)
        
        # Update the shared x_range
        self.shared_x_range.start = start
        self.shared_x_range.end = end
    
    def _figure(self, shared_x_range):
        """Create and configure the figure"""
        fig = bplt.figure(width_policy="max", height=self.height, tools='', x_range=shared_x_range)

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
        # Add multi_line renderer for data
        self.figure.multi_line(
            xs='xs', ys='ys', line_color='colors', line_alpha=0.3, source=self.source_lines
        )

class Pileupy:
    def __init__(self, region, genome=None, reference=None, control=False):
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
        
        self.control = control
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
            self.genome = None
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
        self._create_base_tracks()
        # if genome:
        # else:
        #     self._create_base_tracks()

    def _create_base_tracks(self):
        """Create the basic genome browser tracks"""
        print('here')
        if self.genome:
            print('is genome')
            self.ref = pysam.FastaFile(self.genome_paths[self.genome]['ref'])
            self.ref_seq = self.ref.fetch(self.chrom, self.start, self.end)
            self.genes = load_refgene(self.genome_paths[self.genome]['refgene'])
            self.cytobands, self.chr_lengths = load_cytoband(self.genome_paths[self.genome]['cytoband'])
        else:
            print('is reference')
            self.ref = pysam.FastaFile(self.reference)
            self.genes = []
            self.cytobands = []
            self.chr_lengths = {}
        
        self.shared_x_range = Range1d(self.start, self.end)

        # self.track_control()
        # self.track_idiogram()
        if self.genome:
            if self.control:
                self.add_track_control()
            self.add_track_idiogram()
            self.add_track_reference()
            self.add_track_gene()
        else:
            if self.control:
                self.add_track_control()
            self.add_track_reference()
        
        print(self.tracks)
        print('track made')

        # self._update_layout()
        # self._update_tracks(self.chrom, self.start, self.end)

    def add_track_control(self):
        """Create control track"""
        track_control = TrackControl(self.chr_lengths, self.chrom, self.start, self.end, self._update_tracks)
        self.tracks.append(track_control)

    def add_track_idiogram(self):
        """Create chromosome idiogram track with cytobands"""
        track_idiogram = TrackIdiogram(self.cytobands, self.chr_lengths, self.chrom, self.start, self.end, self.shared_x_range)
        self.tracks.append(track_idiogram)
    
    def add_track_reference(self):
        """Create the reference sequence track"""
        # Create the figure first with shared x range
        track_reference = TrackReference(self.ref, self.chrom, self.start, self.end, self.shared_x_range)
        self.tracks.append(track_reference)
    
    def add_track_gene(self):
        """Add an annotation track from BED file
        
        Args:
            bed_path: Path to BED file
        """
        self.genes
        track_gene = TrackGene(self.genes, self.chrom, self.start, self.end, self.shared_x_range)
        self.tracks.append(track_gene)

    def add_track_alignment(self, bam_path, mode='collapsed', height=400):
        """Add an alignment track to the visualization
        
        Args:
            bam_path: Path to BAM file
            mode: Display mode ('collapsed' or 'expanded')
        """
        bam = pysam.AlignmentFile(bam_path)
        track_alignment = TrackAlignment(self.ref, bam, self.chrom, self.start, self.end, self.shared_x_range, height=height, update_callback=self._update_tracks)
        self.tracks.append(track_alignment)

    def add_track_annotation(self, bed_path):
        """Add an annotation track from BED file
        
        Args:
            bed_path: Path to BED file
        """
        track_annotation = TrackAnnotation(bed_path, self.chrom, self.start, self.end, self.shared_x_range)
        self.tracks.append(track_annotation)

    def add_track_df(self, df, height=100):
        """Add an annotation track from BED file
        
        Args:
            bed_path: Path to BED file
        """
        track_df = TrackDataFrame(df, self.chrom, self.start, self.end, self.shared_x_range, height=height)
        self.tracks.append(track_df)

    def _update_tracks(self, chrom, start, end, filter=[]):
        """Update the tracks after adding new tracks"""
        # Update the shared x_range
        self.shared_x_range.start = start
        self.shared_x_range.end = end
        self.filter = filter
        
        if self.control:
            for track in self.tracks[1:]:
                track._update(self.ref, chrom, start, end, self.shared_x_range)
        else:
            for track in self.tracks:
                track._update(self.ref, chrom, start, end, self.shared_x_range)
        
        print('updated')

    def show(self):
        """Display the visualization"""
        bplt.show(column(
                [track.figure for track in self.tracks],
                sizing_mode="stretch_both",  # Make layout responsive in both directions
                width_policy="max",
                height_policy="max"
            ))
    
    def serve(self, port=5006, num_procs=1, show=True):
        """
        Serve the Pileupy browser as a Bokeh application.
        
        :param port: Port to serve the application on.
        :param num_procs: Number of worker processes.
        :param show: Whether to open the browser automatically.
        """
        def bkapp(doc):
            # Create a fresh layout for each document
            print('serving')
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