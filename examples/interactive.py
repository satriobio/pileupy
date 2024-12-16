from pileupy.main import Pileupy

browser = Pileupy('chr22:24376166-24376456', genome='hg19')
# chr22:24076456-24376456
browser.track_alignment('data/gstt1_sample.bam')
browser.track_annotation('data/mod.bed')

browser.serve()