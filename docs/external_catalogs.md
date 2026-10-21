# External catalogs (F2.1 / F2.2 / F2.3)

The novelty claim needs at least one pinned external lncRNA catalog. Without one,
the pipeline still runs and records `novelty_status = not_evaluated`.

## What is needed

One **BED per source** (BED6 or BED12), same assembly as the reference
(`hg38`/`GRCh38`) and the same chromosome naming as the reference GTF
(`1`, `2`, … as in the Ensembl/GENCODE GTF used by the run; **no `chr` prefix**
if the GTF has none). Column 1 is the source feature id.

| name | source | required |
|---|---|---|
| `lncipedia` | LNCipedia (hg38) | yes |
| `noncode` | NONCODE (hg38) | yes |
| `rnacentral` | RNAcentral | optional |
| `fantomcat` | FANTOM CAT | optional |

## How to fill the spec

Edit `docs/external_catalogs.json` and set `version`, `url` and (once downloaded)
`file` for each source. Record the exact release — the point of the spec is that
someone else can reproduce the overlap.

```json
{"name": "lncipedia", "version": "<release>", "url": "<where you got it>",
 "file": "assets/catalogs/lncipedia.hg38.bed"}
```

Download from the project sites (the researcher picks and pins the release):

- LNCipedia: <https://lncipedia.org> — download the hg38 transcript set and
  convert the GTF/FASTA coordinate table to BED.
- NONCODE: <http://www.noncode.org> — download the human hg38 BED/GTF and convert
  to BED if needed.

Conversion to BED6 keeps it simple:

```bash
# from a GTF
awk -F'\t' '$3=="transcript" {
  match($9, /transcript_id "[^"]+"/); id=substr($9, RSTART+15, RLENGTH-16);
  print $1"\t"$4-1"\t"$5"\t"id"\t0\t"$7 }' catalog.gtf > catalog.bed
```

## Register and check

`bin/external_catalogs.py` has two modes; `check` never downloads.

```bash
python3 bin/external_catalogs.py --spec docs/external_catalogs.json --mode check \
    --out run_info/external_catalogs.tsv
```

Then pass the local BEDs to the run:

```
--lncipedia_bed assets/catalogs/lncipedia.hg38.bed \
--noncode_bed   assets/catalogs/noncode.hg38.bed
```

`run_info/external_catalogs.tsv` records feature counts and md5. The F2.3 guard
(`CHECK_NOVELTY`) then fails fast if a BED has 0 features, shares no chromosome
with the consensus annotation, or yields a 0 hit rate; otherwise it records
`novelty_status = evaluated`.
