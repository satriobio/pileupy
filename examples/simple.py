from pileupy.main import Pileupy

browser = Pileupy('chr22:24376166-24376456', genome='hg19')
browser.add_track_alignment('data/gstt1_sample.bam')

browser.show()
