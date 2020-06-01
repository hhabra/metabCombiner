##This file contains utility functions for matching character string values

#' Ignore Bracketed Strings
#'
#' @description
#' This function is generally used to ignore certain strings that are
#' bracketed for the purposes of string comparisons
#'
#' @param s  a character vector
#'
#' @param brackets  any one or a combination of bracket characters
#'
#' @return character vector with bracketed strings replaced with empty string
#'
#' @noRd
ignore_brackets = function(s, brackets)
{
    if("(" %in% brackets | ")" %in% brackets | "()" %in% brackets)
        s = gsub("^\\([^&]*\\)$", "", s)

    if("[" %in% brackets | "]" %in% brackets | "[]" %in% brackets)
        s = gsub("^\\[[^&]*\\]$", "", s)

    if("{" %in% brackets | "}" %in% brackets | "{}" %in% brackets)
        s = gsub("^\\{[^&]*\\}$", "", s)

    return(s)
}


#' Determine Matching Identity Strings
#'
#' @param s1  One of two character vectors to be compared
#'
#' @param s2  One of two character vectors to be compared
#'
#' @param match The value to give for matching s1 and s2 strings
#'
#' @param mismatch The value to give for mismatching or empty s1 and s2 strings
#'
#' @param brackets  bracketed character strings of these types will be
#' ignored according to this argument
#'
#' @return  vector. Indices of matching strings are equal to value of match,
#'          with mismatching strings equal to value of mismatch
#'
#' @noRd
compare_strings <- function(s1,s2, match, mismatch, brackets)
{
    s1 = ifelse(is.na(s1), "", s1)
    s2 = ifelse(is.na(s2), "", s2)

    s1 = ignore_brackets(s1, brackets)
    s2 = ignore_brackets(s2, brackets)

    ifelse(!(s1 == "" | s2 == "") & tolower(s1) == tolower(s2), match, mismatch)
}

#' Determine Matching Adduct Strings
#'
#' @param s1  One of two character vectors to be compared
#'
#' @param s2  One of two character vectors to be compared
#'
#' @param match The value to give for matching s1 and s2 strings
#'
#' @param mismatch The value to give for mismatching or empty s1 and s2 strings
#'
#' @param brackets  bracketed character strings of these types will be
#' ignored according to this argument
#'
#' @return  character vector. Indices of mismatching strings are equal to value
#'          of mismatch, with remaining indices left as empty characters.
#'
#' @noRd
compare_adduct_strings <- function(s1,s2, match, mismatch, brackets)
{
    s1 = ifelse(is.na(s1), "", s1)
    s2 = ifelse(is.na(s2), "", s2)

    s1 = ignore_brackets(s1, brackets)
    s2 = ignore_brackets(s2, brackets)

    ifelse(!(s1 == "" | s2 == "") & tolower(s1) != tolower(s2), mismatch, match)
}


