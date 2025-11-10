# =============================================================================
# DATASETS FOR MCC PACKAGE
# =============================================================================

#' Titanic Dataset
#'
#' Classical Titanic passenger data containing information on survival
#' status according to Class and Sex. Frequencies have been expanded so that
#' each row represents one passenger.
#'
#' @format A data frame with columns:
#' \describe{
#'   \item{Class}{Passenger class (1st, 2nd, 3rd, Crew)}
#'   \item{Sex}{Gender (Male, Female)}
#'   \item{Survived}{Survival status (Yes, No)}
#' }
#'
#' @examples
#' data(titanic)
#' head(titanic)
#'
"titanic"


#' UCBAdmissions dataset (Berkeley)
#'
#' Admissions information for the six largest departments of the University
#' of California, Berkeley, in Fall 1973. Frequencies have been expanded so that
#' each row represents one applicant.
#'
#' @format A data frame with columns:
#' \describe{
#'   \item{Dept}{Department (A–F)}
#'   \item{Gender}{Applicant gender (Male, Female)}
#'   \item{Admit}{Admission status (Admitted, Rejected)}
#' }
#' @source \url{https://stat.ethz.ch/R-manual/R-devel/library/datasets/html/UCBAdmissions.html}
#' @examples
#' data(berkley)
#' head(berkley)
#'
"berkley"
