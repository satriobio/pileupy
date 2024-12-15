# PileuPy

A Python package to build interactive, rich pileup easily. just add whatever track(s) you need! 

## Quick start

Human genome 
```
# Initialize browser
browser = Pileupy('chr22:24376166-24376456', genome='hg19')

# Add tracks
browser.track_alignment('gstt1_sample.bam')
browser.track_annotation('mod.bed')
```

Custom reference 
```
# Initialize browser
browser = Pileupy('chr22:24376166-24376456', reference='ref.fa')

# Add tracks
browser.track_alignment('gstt1_sample.bam')
browser.track_annotation('mod.bed')
```

Web server
```
pileupy serve
```

## Installation

## Interactive view

## Creating report

## Examples

## License


sv's 
http://localhost:8080/examples/structural-variants.html

quick-start

URL: https://www.encodeproject.org/files/ENCFF716VWO/@@download/ENCFF716VWO.bigWig
URL: https://www.encodeproject.org/files/ENCFF669DTI/@@download/ENCFF669DTI.bigWig
URL: https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/ncbiRefSeqSelect.txt.gz

vcf and cram
https://igv.org/web/release/2.15.5/examples/events.html

phasing


copy number
https://igv.org/web/release/2.15.5/examples/copyNumber.html

events
https://igv.org/web/release/2.15.5/examples/events.html

splice junction
Mutation Annotation Format
https://igv.org/web/release/2.15.5/examples/maf-tcga.html

bigInteract
https://igv.org/web/release/2.15.5/examples/interact.html

loading igv sessions and files from html links

coloring GFF and GTF annotation tracks
https://igv.org/web/release/2.15.5/examples/gff-colors.html

