# PangeBlocks-wide
PangeBlocks for Human Genomes

### How to Use

| Flag            | Description                   |
| --------------- | ----------------------------- |
| `--fasta`, `-f` | Path to the input FASTA file  |
| `-k`            | k-mer size for Jaccard index  |

The sequences are expected to be given in a single FASTA file like this:
```txt
>seq1
ATCGATCGATCG
>seq2
GCTAGCTA
>seq3
TTATCGATCGAATC
```

The graph is produced on stdout with the gfa format.
