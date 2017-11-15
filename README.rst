PPanGGOLiN : Depicting microbial species diversity via a partitioned pangenome graph
========================================================

This tools compile gene content of a taxonomic unit (pangenome) and partition it into 3 partitions (persistent, shell, cloud) via a Bernoulli Mixture Model optionally smoothed using a Markov Random Field.

Moreover the tool build a graph of gene neighborhood and project the partion on it to obtain the partitioned pangenome graph.

It can be instaled via pip `pip install ppanggolin`

#Input


Example of usage:
```
>>> from ppanggolin import PPanGGOLiN
>>> PPanGGOLiN()
```
