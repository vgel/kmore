# kmore

kmer matching experiments.

build a match table:

```bash
cargo run --release build-wtab data/wtab/ data/wtab_18mer.dat
```

this will always generate a fixed-size match table (8.1g), with varying fill rates. low-fill-rate tables compress well (but must be decompressed to use).

use a match table:

```bash
$ cargo run --release -- query data/wtab/GCF_000864105.1_Influenza_A.fna data/wtab_18mer.dat
(start querying...)
(querying took 0.01s)
matched 26524 / 26524 windows (100%)

$ cargo run --release -- query data/smalltomato.fna data/wtab_18mer.dat
(start querying...)
(querying took 0.05s)
matched 0 / 137834 windows (0%)
```


