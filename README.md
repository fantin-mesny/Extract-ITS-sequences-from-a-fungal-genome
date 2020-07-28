# Extract-ITS-sequences-from-a-fungal-genome

ITSx (Bengtsson‐Palme *et al.*, 2013) is still the reference method to extract ITS sequences from genomic fasta files. It is however really slow.
More recently, Barrnap - a fast and accurate method to identify the location of ribosomal RNA genes - was developed.

By combining these two softwares and performing sequences comparison, this script allows the fast extraction of ITS sequences from fungal genomes.

```
python extractITS.py -which ITS2 -i genome.fasta -o ./output/ -name mySpecies
```

## Dependencies

- Python 3.x
- Barrnap: https://github.com/tseemann/barrnap
- ITSx: https://microbiology.se/software/itsx
- Biopython
- Pandas
