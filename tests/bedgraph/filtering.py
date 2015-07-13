import pybio
import sys

bed = pybio.data.Bedgraph()
bed.load("filtering.bed", id="a", meta=["filtering.bed"])
bed.filter()
bed.save("filtering.res.bed")
