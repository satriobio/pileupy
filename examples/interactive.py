from pileupy.main import Pileupy

# browser = Pileupy('chr22:24376166-24376196', genome='hg19', control=True)
browser = Pileupy('chr22:24376166-24376196', reference='data/hg19.fa', control=True)
browser.add_track_alignment('data/gstt1_sample.bam')

browser.serve()