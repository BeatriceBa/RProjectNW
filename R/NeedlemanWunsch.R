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
    matrix <- matrix(0, nrow = nchar(sequenceA)+1, ncol = nchar(sequenceB)+1)

    for(i in 2:nrow(matrix)) {
        matrix[i,1] <- gapPenalty * (i-1)
    }

    for(j in 2:ncol(matrix)) {
        matrix[1,j] <- gapPenalty * (j-1)
    }

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
    matrixDirection <- matrix
    alignmentA <- ""
    alignmentB <- ""
    i <- nrow(matrixDirection)
    j <- ncol(matrixDirection)
    while (i>1 | j>1) {
        matrixDirection[i,j] <- "*"

        if(isItAMatch(sequenceA,sequenceB,i,j)) {
            alignmentA <- paste(substring(sequenceA, i-1, i-1),alignmentA)
            alignmentB <- paste(substring(sequenceB, j-1, j-1),alignmentB)
            i <- i-1
            j <- j-1
        }
        else {
            top <- matrix[i-1,j]
            left <- matrix[i,j-1]
            diagonal <- matrix[i-1,j-1]

            if (top == max(c(top,left,diagonal))) {
                alignmentA <- paste(substring(sequenceA, i-1, i-1),alignmentA)
                alignmentB <- paste("-",alignmentB)
                i <- i-1
            }

            else if (left == max(c(top,left,diagonal))){
                alignmentA <- paste("−",alignmentA)
                alignmentB <- paste(substring(sequenceB, j-1, j-1),alignmentB)
                j <- j-1
            }

            else  {
                alignmentA <- paste(substring(sequenceA, i-1, i-1),alignmentA)
                alignmentB <- paste(substring(sequenceB, j-1, j-1),alignmentB)
                i <- i-1
                j <- j-1
            }
        }
    }
    result <- list(matrixDirection = matrixDirection,
                   alignmentA = alignmentA,
                   alignmentB = alignmentB)
    return(result)
}

#' Needleman Wunsch
#'
#' @description
#' Core function of the package, performs the Needleman-Wunsch algorithm
#'
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
    scoresMatrix <- creatingMatrix(matchCost, mismatchCost, gapPenalty, sequenceA, sequenceB)
    alignment <- performTraceback(scoresMatrix,sequenceA,sequenceB)

    result <- list(scoresMatrix = scoresMatrix,
                   matrixDirection = alignment$matrixDirection,
                   alignmentA = alignment$alignmentA,
                   alignmentB = alignment$alignmentB)
    return(result)
}
