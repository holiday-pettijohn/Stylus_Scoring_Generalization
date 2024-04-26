# Introduction to Stylus

## What is Stylus?

Stylus is a virtual world model for genetic inheritance and mutation using Chinese characters as representationsn of genomes, as drawn compared to their ideal 'archetypes', which are canonical Chinese characters composed of a certain set of strokes. More legible characters are considered better 'fit' genomes, whereas less readable ones are considered less fit on a gradual scale that ranges from slight imperfection to illegibility.

## Stroke Matching


The algorithm which measures the fitness of a specific gene by comparing it to its respective archetype by using the stroke data by which both the gene and archetype are encoded. However, this algorithm can only score genes against archetypes when it has an appropriate mapping of which strokes should be paired with which genes. However, there are O(n!) possible mappings between any n archetype and gene strokes, and the canonical scoring algorithm becomes infeasibly expensive at higher stroke counts. Our task is to use a heuristic scoring algorithm to search for the optimal gene stroke to archetype pairing, that is, the one which returns the highest canonical score, in an acceptable amount of time (i.e. not O(n!), or ideally anything non-polynomial).

## Matching Algorithm Tutorial

The current greedy matching algorithm matches strokes to their archetype counterpart by approximating the similarity of each of the gene and archetype strokes, then pairing them together based on the one with the best score in the order of archetype stroke length.

Currently, error between strokes is defined as the maximum deviation between a stroke and its archetype counterpart at any of its critical points (that is, the segment edges at which the stroke's geometry is defined).

This algorithm is invoked in the Heuristic_Exepriment notebook.
