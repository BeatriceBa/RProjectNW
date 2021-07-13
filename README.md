# Needleman Wunsch algorithm 

## Description

Sequence alignment is a way of arranging sequences of DNA, RNA, or protein to identify regions of similarity that may be a consequence of functional, structural, or evolutionary relationships between the sequences. The alignment of exactly two sequences is called pairwise alignment, and it can be performed globally or locally. The Needleman-Wunsch package provides functions to perform a pairwise global sequence alignment using the Needleman-Wunsch algorithm.

## Functions

### Align two sequences

Aligning two sequences using the Needleman-Wunsch algorithm is now an easy task. With the globalAlignnmentNeedlemanWunsch function. What the user needs are the costs of matches and mismatches, the gap penalty and the two sequences.
It will return a list containing the scores matrix, the direction matrix and the alignment

`print(globalAlignnmentNeedlemanWunsch(7,-3,-4,"GTT","GCATT"))`

### Creating the scores matrix

The first step of the Needleman-Wunsch algorithm is the creation of scores matrix. The function creatingMatrix computes this matrix. As for the globalAlignmentNeedlemanWunsch function, what the user needs are the costs of matches and mismatches, the gap penalty and the two sequences.

`print(creatingMatrix(7,-3,-4,"GTT","GCATT"))`

### Traceback

In order to align two sequences, the creation of the scores matrix is followed by the traceback step. The function performTraceback generates a matrix with the directions of the alignment while computing the alignment itself. It takes as input a scores matrix and two sequences.

`performTraceback(scoresMatrix, "GTT", "GCATT")`
