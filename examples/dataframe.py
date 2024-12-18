from pileupy.main import Pileupy
import pandas as pd
import numpy as np

# Base position input
contig = "chr22"
start = 24376166
end = 24376456

# Generate positions
positions = np.arange(start, end + 1)

# Number of reads you want to generate
num_reads = 5  # You can change this number based on the number of reads you want

# Create an empty list to store data for each read
data = []

# Generate Gaussian noise for each read
for read_id in range(num_reads):
    values = np.random.normal(loc=0, scale=1, size=len(positions))  # Gaussian noise for each read
    read_data = {
        'chrom': [contig] * len(positions),
        'start': positions,
        'end': positions + 1,
        'value': values,
        'read_names': [read_id] * len(positions)
    }
    df_read = pd.DataFrame(read_data)
    data.append(df_read)

# Concatenate all reads' data into a single DataFrame
final_df = pd.concat(data, ignore_index=True)

browser = Pileupy('chr22:24376166-24376456', genome='hg19')
browser.add_track_alignment('data/gstt1_sample.bam')
browser.add_track_df(final_df)

browser.serve()