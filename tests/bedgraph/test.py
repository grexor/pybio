import pybio
import sys

bed = pybio.data.Bedgraph()
bed.load("file_a.bed", meta=["file_a.bed"])
bed.load("file_b.bed", meta=["file_b.bed"])
bed.cluster()

bed.save("file_c.bed", db_filter="present", db_save="raw", min_cDNA=1)
bed.save("file_c.scale.bed", db_filter="raw", db_save="scale", min_cDNA=1)
bed.save("file_c.present.bed", db_filter="raw", db_save="present", min_cDNA=1)
bed.save("file_c.meta.bed", db_filter="raw", db_save="meta", min_cDNA=1)
