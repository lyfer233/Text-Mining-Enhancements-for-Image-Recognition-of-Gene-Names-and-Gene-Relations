# Gene Relation Enhancement Results

### Segmentation

First we need to complete the three levels for an article, see the [cutting_full_text folder](./cutting_full_text) for the implementation of this piece.

### Generate gene occurrence

We calculate the ground truth and the gene occurrence of each level. [For more details](./calc_relation/gene_relation_comparison.py)

### Calculate the confusion matrix

We can use this [python file](./calc_relation/gene_relation_analysis.py) to calculate the obfuscation matrix to calculate the precision, recall and F1