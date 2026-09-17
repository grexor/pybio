# CLIP & interval data

`pybio.data.Bedgraph` is a class for working with per-position signal data keyed by chromosome, strand and position — the kind of sparse, strand-specific counts produced by CLIP-seq crosslink sites, 3'-end sequencing, or any other bedGraph-style track. It loads bedGraph files into an in-memory position → value map and gives you clustering, region queries and normalization on top.

## Loading data

```python
import pybio.data
bg = pybio.data.Bedgraph("sample.bed")
```

`Bedgraph(filename)` loads immediately; you can also construct an empty one and call `.load(filename)` one or more times — repeated loads accumulate onto the same object, which is useful for merging replicates:

```python
bg = pybio.data.Bedgraph()
bg.load("replicate1.bed")
bg.load("replicate2.bed")  # counts are added, not replaced
```

Each line of the input is expected in bedGraph format (`chr  start  stop  value`, optionally gzipped, `track` and `#` lines ignored). The sign of `value` determines strand (`+` for ≥0, `-` for negative) unless you pass `force_strand`. `min_cDNA` drops positions below a minimum raw count while loading.

## Querying values

```python
bg.get_value("1", "+", 11012344)               # value at a single position
bg.get_region("1", "+", 11012344, -20, 20)      # summed value over a window around a position
bg.get_vector("1", "+", 11012300, 11012400)     # per-position values across a range, as a list
```

`get_vector()` reverses the returned list on the minus strand, so it always reads 5'→3' along the transcript — the same convention as [`seq()`](sequences.md#retrieving-genomic-sequence).

To iterate over every stored position:

```python
for chrom, strand, pos, value in bg.fetch():
    ...
```

## Clustering positions

Real crosslink/end-sequencing data is noisy — nearby positions on the same strand often represent the same underlying site. `cluster()` collapses that: for each chromosome/strand, it walks positions from highest value to lowest, and for each one sums every position within `region_up` upstream / `region_down` downstream into it, removing the positions it consumed:

```python
bg.cluster(region_up=150, region_down=150)
```

`filter(min_distance=125)` is a lighter alternative: instead of summing neighbors in, it just keeps the highest-value position in every `min_distance`-nt window and discards the rest.

## Normalizing and saving

```python
bg.norm()  # computes CPM (count per million) into bg.cpm, based on bg.total_raw
bg.save("clustered.bed", db_save="raw", min_raw=2)
```

`save()` writes the object back out in bedGraph format; `db_save` picks which internal table (`"raw"` or `"cpm"`) to write, and `min_raw`/`min_cpm`/`min_support` filter which positions are included.

## Converting to BigWig

Once you have a final bedGraph file, `pybio.data.bedgraph_bigwig(filename_bed, filename_bw, genome)` wraps `bedGraphToBigWig` to produce a BigWig track (requires `bedGraphToBigWig` on your `PATH`, and a `.chrs` chromosome-sizes file for `genome`). Because BigWig can't represent both strands in one track, split plus/minus strand data into separate bedGraph files before converting.
