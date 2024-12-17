# Examples

## Basic setup

This setup will generate standalone HTML report and open new tab in browser.

```
from pileupy.main import Pileupy

browser = Pileupy('chr22:24376166-24376456', genome='hg19')
browser.add_track_alignment('data/demo.bam')
browser.show()
```

This setup will start interactive genome browser.

```
from pileupy.main import Pileupy

browser = Pileupy('chr22:24376166-24376456', genome='hg19', control=True)
browser.add_track_alignment('data/demo.bam')
browser.serve()
```

You can also start interactive genome browser with command-line.
```
pileupy view
```

## Tracks

### Reference

Choose available genome reference detailed bellow using `genome` or specify reference path using `reference`. Choosing reference from available genome will also load Idiogram and RefSeq track. 

- Human hg19 `hg19`
- Human hg38 `hg38`

```
from pileupy.main import Pileupy

browser = Pileupy('chr22:24376166-24376456', genome='hg19')
# browser = Pileupy('chr22:24376166-24376456', reference='/foo/bar/hg19.fa')
```

### Alignment

Add aligment track using `add_track_alignment`. Format supported is `SAM`, `BAM`, and `CRAM`.

```
from pileupy.main import Pileupy

browser = Pileupy('chr22:24376166-24376456', genome='hg19')
browser.add_track_alignment('data/demo.bam')
```

Disable SNP highlight

```
browser.add_track_alignment('data/demo.bam', snp=False)
```

Group reads by phasing

```
browser.add_track_alignment('data/demo.bam', phasing=True)
```

### Annotations

Add annotation track using `add_track_annotation`. Format supported is `BED`.

```
from pileupy.main import Pileupy

browser = Pileupy('chr22:24376166-24376456', genome='hg19')
browser.add_track_alignment('data/demo.bam')
```

### Data Frame

Add data frame track using `add_track_annotation`. Input used is pandas data frame with at least following columns `chrom`, `start`, `value`

```
from pileupy.main import Pileupy

browser = Pileupy('chr22:24376166-24376456', genome='hg19')
browser.add_track_alignment('data/demo.bam')
```