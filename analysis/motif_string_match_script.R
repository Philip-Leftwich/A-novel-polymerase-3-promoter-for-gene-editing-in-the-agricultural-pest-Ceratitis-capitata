library(tidyverse)
library(stringdist)


seq_7sk <- read_csv("data/raw/seqs_7sk.csv")

seq_7sk <- seq_7sk |>
  mutate(
    `7sk_rna` = str_replace_all(`7sk_rna`, " ", ""),
    `7sk_promoter` = str_replace_all(`7sk_promoter`, " ", "")
  )

last_match_start <- function(string, pattern) {
  map_int(
    str_locate_all(string, pattern),
    ~ {
      if (nrow(.x) == 0) NA_integer_ else .x[nrow(.x), "start"]
    }
  )
}

last_match_string <- function(string, pattern) {
  map_chr(
    str_extract_all(string, pattern),
    ~ {
      if (length(.x) == 0) NA_character_ else .x[length(.x)]
    }
  )
}

# For each string, if no regex match, find closest substring to consensus
closest_partial <- function(string, consensus) {
  motif_len <- nchar(consensus)
  map_chr(
    string,
    ~ {
      if (nchar(.x) < motif_len) {
        return(NA_character_)
      }
      # Extract all substrings of motif length
      starts <- seq_len(nchar(.x) - motif_len + 1)
      substrings <- substring(.x, starts, starts + motif_len - 1)
      # Return substring with minimum edit distance to consensus
      dists <- stringdist(substrings, consensus, method = "hamming")
      substrings[which.min(dists)]
    }
  )
}

# ...existing code...

# Returns start position of closest substring match within window
closest_partial_start <- function(string, consensus) {
  motif_len <- nchar(consensus)
  map_int(
    string,
    ~ {
      if (nchar(.x) < motif_len) {
        return(NA_integer_)
      }
      starts <- seq_len(nchar(.x) - motif_len + 1)
      substrings <- substring(.x, starts, starts + motif_len - 1)
      dists <- stringdist(substrings, consensus, method = "hamming")
      starts[which.min(dists)]
    }
  )
}

string_match <- function(seq_7sk) {
  PSEA_pattern <- "TAATTCCCAAGTGCTTATTTG"
  TATA_pattern <- "(A|C|T)(A|T)TA(A|T)A"

  seq_7sk |>
    mutate(
      seq_len = nchar(`7sk_promoter`),
      tata_window = substr(`7sk_promoter`, seq_len - 50 + 1, seq_len - 20 + 1),
      psea_window = substr(`7sk_promoter`, seq_len - 100 + 1, seq_len - 30 + 1),
      TATA_win_pos = last_match_start(tata_window, TATA_pattern),
      PSEA_win_pos = last_match_start(psea_window, PSEA_pattern),
      # Record match type before filling NAs
      PSEA_match_type = if_else(is.na(PSEA_win_pos), "partial", "exact"),
      TATA_match_type = if_else(is.na(TATA_win_pos), "partial", "exact"),
      # Fill NA positions with fuzzy match position
      TATA_win_pos = if_else(
        is.na(TATA_win_pos),
        closest_partial_start(tata_window, TATA_pattern),
        TATA_win_pos
      ),
      PSEA_win_pos = if_else(
        is.na(PSEA_win_pos),
        closest_partial_start(psea_window, PSEA_pattern),
        PSEA_win_pos
      ),
      TATA = if_else(
        !is.na(TATA_win_pos),
        50L - TATA_win_pos + 1L,
        NA_integer_
      ),
      Insect_PSEA = if_else(
        !is.na(PSEA_win_pos),
        100L - PSEA_win_pos + 1L,
        NA_integer_
      ),
      TATA_match = last_match_string(tata_window, TATA_pattern),
      PSEA_match = last_match_string(psea_window, PSEA_pattern),
      TATA_match = if_else(
        is.na(TATA_match),
        closest_partial(tata_window, TATA_pattern),
        TATA_match
      ),
      PSEA_match = if_else(
        is.na(PSEA_match),
        closest_partial(psea_window, PSEA_pattern),
        PSEA_match
      )
    ) |>
    select(
      species,
      seq_id,
      seq_len,
      Insect_PSEA,
      PSEA_match,
      PSEA_match_type,
      TATA,
      TATA_match,
      TATA_match_type
    )
}
string_match(seq_7sk)
