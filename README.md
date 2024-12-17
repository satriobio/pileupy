![logo](docs/img/icon.png)

# PileuPy

A python package for visualizing genomic data easily. 

## Quick start

Generate standalone HTML report
```
from pileupy.main import Pileupy

browser = Pileupy('chr22:24376166-24376456', genome='hg19')
browser.add_track_alignment('gstt1_sample.bam')
browser.show()
```

Start interactive browser
```
from pileupy.main import Pileupy

browser.add_track_alignment('gstt1_sample.bam')
browser.add_track_annotation('mod.bed')
browser.serve()
```

Start interactive browser using command-line

```
pileupy --region chr22:24376166-24376456 --genome hg19  --alignment gstt1_sample.bam
```

## Installation

```
pip install pileupy
```

## License

GNU General Public License v3.0