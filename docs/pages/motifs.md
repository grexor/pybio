# Sequence & motif tools

`pybio.sequence` is a small toolkit for working with raw sequence strings: reverse complementation, IUPAC-aware motif search, signal smoothing, and nucleotide composition plots. It's independent of genome downloads — everything here works on plain Python strings.

## Reverse complement

```python
import pybio
pybio.sequence.reverse_complement("ATCGN")
# 'NCGAT'
```

Handles the full IUPAC ambiguity alphabet (`R`, `Y`, `K`, `M`, `S`, `W`, `B`, `D`, `H`, `V`, `N`) and preserves letter case.

## Motif search

`search(input_string, motif_list)` scans a sequence for one or more motifs and returns both the list of match spans and a 0/1 coverage vector the length of the input:

```python
matches, coverage = pybio.sequence.search("AAATGCATTTGCA", ["ATG", "TGC"])
```

`motif_list` can be a Python list of motifs, or a single string with motifs separated by `_` (e.g. `"ATG_TGC"`). Motifs whose first and last character are one of the IUPAC ambiguity codes `R`, `Y`, `S`, `W` are automatically expanded to every matching A/T/C/G combination before searching — pass `expand()` a motif directly if you want to see the expansion:

```python
pybio.sequence.expand("RAAY")
# ['AAAC', 'AAAT', 'GAAC', 'GAAT']
```

Pass `strict=True` to `search()` to require every motif in `motif_list` to match at least once, returning no matches at all otherwise.

## Smoothing a signal vector

Three helpers operate on a numeric vector (typically the coverage vector from `search()`, or any other per-position signal):

- `filter(vector, hw=30, hwt=1)` — convolves with a `2*hw+1`-wide flat window, then thresholds back to 0/1 at `hwt`. Useful for turning sparse motif hits into contiguous "hit regions" within `hw` bases of each other.
- `convolve(vector, hw)` — the same convolution, without thresholding — smooths a signal without binarizing it.
- `extend(vector, window_size)` — sums each position over a `window_size`-wide neighborhood, without normalizing.

## Nucleotide composition plots

`draw(sequences, fname)` plots the per-position A/U/C/G percentage across a list of same-purpose sequences (e.g. sequences centered on a motif or crosslink site) and saves it to `fname`:

```python
pybio.sequence.draw(["ATGCATGC", "ATGGATCC", "ATGCATCC"], "composition.png")
```

Sequences of unequal length are padded to the longest one; positions past the end of a shorter sequence are excluded from that sequence's contribution.
