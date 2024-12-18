![logo](img/icon.png)

# PileuPy

A python package for visualizing genome browser tracks easily. Built on top of Bokeh to deliver high-performance, customizable visualizations for genomic data.

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

<!-- Start interactive browser using command-line

```
pileupy --region chr22:24376166-24376456 --genome hg19  --alignment gstt1_sample.bam
``` -->
Open the app [http://localhost:5006/](http://localhost:5006/)

![logo](img/interactive.png)

## Installation

```
pip install pileupy
```

## License

GNU General Public License v3.0
