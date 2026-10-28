# Column types of the entity tables.
#
# SQLite keeps a categorical as its labels, a boolean as 0/1 and every integer as INTEGER,
# so the writers record each column's type in `_column_meta` (table_name, column_name,
# dtype, categories): "categorical" or "ordered_categorical" with the categories as a JSON
# array, "bool", an integer dtype ("int64", "Int64", "uint8", ...), a float dtype, or
# "string". Read back here, so a categorical comes back as a factor with its levels in their
# order (unused ones included), a boolean as logical and an integer as integer, as the
# Python reader restores them for pandas.
#
# Categories keep their own JSON type. A categorical's values are its labels, stored as text
# by cytome before 0.3.9 and in the type of its categories since (INTEGER, REAL or TEXT),
# never as codes. String categories are matched to the stored labels as strings; numeric
# ones by value, then by the 15 significant digits of the text SQLite made of a float
# category; boolean ones as booleans (stored 1/0, or as text). A recorded type is applied only where the stored values agree with it:
# `_column_meta` can be stale, since a later column write does not update it, so a value
# outside the recorded categories is kept as an extra level, with a warning, and a column
# whose values do not fit its recorded type is left as stored.

.restore_column_types <- function(x, table, df) {
  if (!ncol(df) || !.has_table(x, "_column_meta")) return(df)
  meta <- DBI::dbGetQuery(x$con,
    "SELECT column_name, dtype, categories FROM _column_meta WHERE table_name = ?",
    params = list(table))
  for (k in seq_len(nrow(meta))) {
    col <- meta$column_name[k]
    if (!col %in% names(df)) next
    v <- .typed_column(df[[col]], meta$dtype[k], meta$categories[k], paste0(table, ".", col))
    if (!is.null(v)) df[[col]] <- v
  }
  df
}

# The column `v` as `dtype` says, or NULL to leave it as stored.
.typed_column <- function(v, dtype, categories, what) {
  if (!length(dtype) || is.na(dtype)) return(NULL)
  if (dtype %in% c("categorical", "ordered_categorical"))
    return(.as_factor(v, categories, dtype == "ordered_categorical", what))
  if (dtype %in% c("bool", "boolean")) return(.as_logical(v, what))
  if (grepl("^u?int(8|16|32|64)$", dtype, ignore.case = TRUE)) return(.as_integer(v))
  if (grepl("^float(16|32|64)$", dtype, ignore.case = TRUE))
    return(if (is.numeric(v) || all(is.na(v))) as.numeric(v) else NULL)
  if (dtype == "string") return(if (is.character(v) || all(is.na(v))) as.character(v) else NULL)
  NULL
}

# Integers that fit R's 32-bit integer become integer; larger ones stay double.
.as_integer <- function(v) {
  if (all(is.na(v))) return(rep(NA_integer_, length(v)))
  if (!is.numeric(v)) return(NULL)
  w <- v[!is.na(v)]
  if (any(w != round(w))) return(NULL)
  if (all(abs(w) <= .Machine$integer.max)) as.integer(v) else as.numeric(v)
}

.logical_values <- function(v) {
  if (is.logical(v)) return(v)
  if (is.numeric(v)) return(ifelse(is.na(v), NA, v != 0))
  t <- tolower(trimws(as.character(v)))
  out <- rep(NA, length(t))
  out[t %in% c("true", "1")] <- TRUE
  out[t %in% c("false", "0")] <- FALSE
  out[!is.na(t) & !t %in% c("true", "1", "false", "0")] <- NA
  attr(out, "unknown") <- !is.na(t) & !t %in% c("true", "1", "false", "0")
  out
}

.as_logical <- function(v, what) {
  if (all(is.na(v))) return(rep(NA, length(v)))
  if (is.numeric(v) && !all(v[!is.na(v)] %in% c(0, 1))) {
    warning("cytome: ", what, " is recorded as boolean but holds other values; left as stored")
    return(NULL)
  }
  out <- .logical_values(v)
  if (any(attr(out, "unknown"))) {
    warning("cytome: ", what, " is recorded as boolean but holds other values; left as stored")
    return(NULL)
  }
  attr(out, "unknown") <- NULL
  out
}

# Level labels for numeric categories: as.character(), or 17 significant digits when two
# categories would otherwise print the same.
.num_levels <- function(cats) {
  l <- as.character(cats)
  if (anyDuplicated(l)) l <- sprintf("%.17g", cats)
  l
}

.as_factor <- function(v, categories, ordered, what) {
  cats <- if (length(categories) && !is.na(categories))
    tryCatch(jsonlite::fromJSON(categories, simplifyVector = TRUE), error = function(e) NULL)
  if (is.null(cats) || is.list(cats)) {
    if (length(categories) && !is.na(categories))
      warning("cytome: the categories recorded for ", what, " cannot be read; it is left as stored")
    return(NULL)
  }
  if (is.logical(cats)) {
    key <- .logical_values(v)
    attr(key, "unknown") <- NULL
    lev <- as.character(cats)
  } else if (is.numeric(cats)) {
    key <- suppressWarnings(as.numeric(v))
    lev <- .num_levels(cats)
  } else {
    key <- as.character(v)
    lev <- as.character(cats)
  }
  codes <- match(key, cats)
  if (is.numeric(cats)) {
    ## Stores written before cytome 0.3.9 kept a float category as the text SQLite makes of
    ## it, 15 significant digits: match those as the Python reader does.
    again <- is.na(codes) & !is.na(key)
    if (any(again)) codes[again] <- match(key[again], signif(cats, 15))
  }
  unknown <- !is.na(v) & is.na(codes)
  if (any(unknown)) {
    extra <- unique(as.character(v[unknown]))
    warning("cytome: ", sum(unknown), " values of ", what, " are not among its recorded ",
            "categories (", paste(utils::head(extra, 3), collapse = ", "),
            if (length(extra) > 3) ", ...", "); they are kept as extra levels")
    codes[unknown] <- length(lev) + match(as.character(v[unknown]), extra)
    lev <- c(lev, extra)
  }
  if (anyDuplicated(lev)) lev <- make.unique(lev)
  structure(as.integer(codes), levels = lev, class = c(if (ordered) "ordered", "factor"))
}
