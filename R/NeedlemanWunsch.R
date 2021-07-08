#' Creation of the scores matrix
#'
#' @description
#' creatingMatrix creates a matrix containing the scores of the alignment for the NW algorithm
#'
#' @param matchCost The cost of a match
#' @param mismatchCost The cost of a mismatch
#' @param gapPenalty The gap penalty(cost of gap/indel)
#' @param sequenceA The first nucleotide sequence (positioned on the rows)
#' @param sequenceB The second nucleotide sequence (positioned on the columns)
#'
#' @return A matrix containing the scores of the alignment
#'
#' @examples
#' sequenceA <- "GTT"
#' sequenceB <- "GCATT"
#' matchCost <- 7
#' mismatchCost <- -3
#' gapPenalty <- -4
#' scoreMatrix <- creatingMatrix(matchCost,mismatchCost,gapPenalty,sequenceA,sequenceB)
#'
#' @export
creatingMatrix <- function(matchCost, mismatchCost, gapPenalty, sequenceA, sequenceB) {
    if (!is.numeric(matchCost)) stop("matchCost must be a number!")

    if (!is.numeric(mismatchCost)) stop("mismatchCost must be a number!")

    if (!is.numeric(gapPenalty)) stop("gapPenalty must be a number!")

    if (typeof(sequenceA)!="character" | typeof(sequenceB)!="character") stop("sequenceA or sequenceB are not strings, you dummy!")

    # Creates a MxN matrix filled with 0's, where
    # M = number of characters of sequenceA + 1
    # N = number of characters of sequenceB + 1
    matrix <- matrix(0, nrow = nchar(sequenceA)+1, ncol = nchar(sequenceB)+1)

    # Fill the first row.
    # It skips the cell [1,1] because there must be a 0 in that cell.
    # The i-1 is because indexes in R starts from 1 instead of 0.
    for(i in 2:nrow(matrix)) {
        matrix[i,1] <- gapPenalty * (i-1)
    }

    # Fill the first column.
    # It skips the cell [1,1] because there must be a 0 in that cell.
    # The j-1 is because indexes in R starts from 1 instead of 0.
    for(j in 2:ncol(matrix)) {
        matrix[1,j] <- gapPenalty * (j-1)
    }

    # Looping from the [2,2] cell fills out the matrix.
    # For each cell it stores the value for the neighbouring cells (top and left) and checks if the
    # top-left diagonal neigbour is either a match or a mismatch, adding, respectively, either the match or the
    # mismatch cost.
    # The current cell assumes the maximum among those 3 values (top,left, top-left diagonal).
    for(i in 2:nrow(matrix)) {
        for(j in 2:ncol(matrix)){
            top <- matrix[i-1,j] + gapPenalty
            left <- matrix[i,j-1] + gapPenalty
            if (isItAMatch(sequenceA,sequenceB,i,j))    diagonal <- matrix[i-1,j-1] + matchCost
            else    diagonal <- matrix[i-1,j-1] + mismatchCost
            matrix[i,j] <- max(top,left,diagonal)
        }
    }
    return(matrix)
}
#' Checks if the current characters are matching
#'
#' @description
#' isItAMatch is a utility function that checks whether the i-th and j-th characters of the two sequences are matching
#'
#' @param sequenceA The first nucleotide sequence (positioned on the rows)
#' @param sequenceB The second nucleotide sequence (positioned on the columns)
#' @param i The index for sequenceA (rows)
#' @param j The index for sequenceB (columns)
#'
#' @return TRUE if the two characters are matching, FALSE otherwise
#'
#' @examples
#' isItAMatch("GTT","GCATT",1,2)
#'
#' @export
isItAMatch <- function (sequenceA, sequenceB, i, j) {
    if (typeof(sequenceA)!="character" | typeof(sequenceB)!="character") stop("sequenceA or sequenceB are not strings, you dummy!")

    if (!is.numeric(i) | !is.numeric(j)) stop("Are you sure you passed the right indexes?")

    characterA <- substring(sequenceA, i-1, i-1)
    characterB <- substring(sequenceB, j-1, j-1)
    if (characterA == characterB) return(TRUE)
    return(FALSE)
}
#' Traceback function
#'
#' @description
#' performTraceback is a function that creates a matrix with the directions of the alignment while computing the alignment itself
#'
#' @param matrix The score matrix computed with the two sequences
#' @param sequenceA The first nucleotide sequence (positioned on the rows)
#' @param sequenceB The second nucleotide sequence (positioned on the columns)
#'
#' @return a list containing the directionMatrix, the first alignment, the second alignment
#'
#' @examples
#' isItAMatch("GTT","GCATT",1,2)
#'
#' @export
performTraceback <- function(matrix, sequenceA, sequenceB) {
    if (typeof(sequenceA)!="character" | typeof(sequenceB)!="character") stop("sequenceA or sequenceB are not strings, you dummy!")

    if (is.null(matrix)) stop("Matrix is null!")

    if (!is.numeric(matrix)) stop("Matrix is not numeric!")

    matrixDirection <- matrix
    alignmentA <- ''
    alignmentB <- ''
    i <- nrow(matrixDirection)
    j <- ncol(matrixDirection)
    # Starting from the bottom-right cell of the matrix, check for the directions
    # Create a matrix that stores the directions of the algorithm, starting from
    # copying the scores matrix and substituting the correct cell with an asterisk
    while (i>1 | j>1) {
        matrixDirection[i,j] <- "*"

        # If the current cells correspond to a match in the sequence, move
        # diagonally and add the corresponding character to both of the sequences
        if(isItAMatch(sequenceA,sequenceB,i,j)) {
            alignmentA <- paste(substring(sequenceA, i-1, i-1),alignmentA)
            alignmentB <- paste(substring(sequenceB, j-1, j-1),alignmentB)
            i <- i-1
            j <- j-1
        }
        # If it is not a match, store the values for the top, left and diagonal
        # neighbouring cells and check which of the aforementioned cells stores
        # the highest value
        else {
            top <- matrix[i-1,j]
            left <- matrix[i,j-1]
            diagonal <- matrix[i-1,j-1]

            # If the neighbouring top cell stores the maximum value,
            # it means top is the next direction.
            # The alignmentA string stores the corresponding character of the
            # string sequenceA (rows), while there is a gap on the alignmentB
            # (sequenceB, columns).
            if (top == max(c(top,left,diagonal))) {
                alignmentA <- paste(substring(sequenceA, i-1, i-1),alignmentA)
                alignmentB <- paste("-",alignmentB)
                i <- i-1
            }
            # If the neighbouring left cell stores the maximum value,
            # it means left is the next direction.
            # The alignmentB string stores the corresponding character of the
            # string sequenceB (columns), while there is a gap on the alignmentA
            # (sequenceA, rows).
            else if (left == max(c(top,left,diagonal))){
                alignmentB <- paste(substring(sequenceB, j-1, j-1),alignmentB)
                alignmentA <- paste("-", alignmentA)
                j <- j-1
            }
            # If the neighbouring diagonal cell stores the maximum value,
            # it means we need to move diagonally.
            # There is a mismatch, both alignmentA and alignmentB store the
            # corresponding characters of sequenceA and sequenceB, respectively.
            else  {
                alignmentA <- paste(substring(sequenceA, i-1, i-1),alignmentA)
                alignmentB <- paste(substring(sequenceB, j-1, j-1),alignmentB)
                i <- i-1
                j <- j-1
            }
        }

    }
    # The result is a named list containing the direction matrix (just to give to the user
    # a graphical visualization of the directions) and the two alignments.
    result <- list(matrixDirection = matrixDirection, alignmentA = alignmentA, alignmentB = alignmentB)
    return(result)
}
#' Needleman Wunsch
#'
#' @description
#' Core function of the package, performs the Needleman-Wunsch algorithm
#'
#'
#' @references \url{https://en.wikipedia.org/wiki/Needleman-Wunsch_algorithm}
#' @param matchCost The cost of a match
#' @param mismatchCost The cost of a mismatch
#' @param gapPenalty The gap penalty(cost of gap/indel)
#' @param sequenceA The first nucleotide sequence (positioned on the rows)
#' @param sequenceB The second nucleotide sequence (positioned on the columns)
#'
#' @return a list containing the score matrix, the direction matrix, the first alignment, the second alignment
#'
#' @examples
#' sequenceA <- "GTT"
#' sequenceB <- "GCATT"
#' matchCost <- 7
#' mismatchCost <- -3
#' gapPenalty <- -4
#' result <- globalAlignnmentNeedlemanWunsch(matchCost,mismatchCost,gapPenalty,sequenceA,sequenceB)
#'
#' @export
globalAlignnmentNeedlemanWunsch <- function(matchCost, mismatchCost, gapPenalty, sequenceA, sequenceB) {
    # Creation of the score matrix
    scoresMatrix <- creatingMatrix(matchCost, mismatchCost, gapPenalty, sequenceA, sequenceB)
    # Traceback with the previously computed matrix
    alignment <- performTraceback(scoresMatrix,sequenceA,sequenceB)

    # The result is a list containing the score matrix, the direction matrix, the first alignment, the second alignment
    result <- list(scoresMatrix = scoresMatrix, matrixDirection = alignment$matrixDirection, alignmentA = alignment$alignmentA, alignmentB = alignment$alignmentB)
    return(result)
}
