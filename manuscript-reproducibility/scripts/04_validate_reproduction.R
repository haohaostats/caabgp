root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
manifest_path <- file.path(root, "checksums", "expected_sha256.csv")
if (!file.exists(manifest_path)) stop("Missing checksum manifest: ", manifest_path)

sha256_file <- function(path) {
  if (.Platform$OS.type == "windows") {
    cmd <- Sys.which("certutil")
    if (!nzchar(cmd)) stop("certutil is required for checksum validation on Windows.")
    out <- system2(cmd, c("-hashfile", shQuote(path), "SHA256"), stdout = TRUE, stderr = TRUE)
    out <- iconv(out, from = "", to = "UTF-8", sub = "")
    hit <- grep("^[0-9A-Fa-f ]{64,}$", trimws(out), value = TRUE)
    if (!length(hit)) stop("Could not compute SHA256 for ", path)
    return(tolower(gsub(" ", "", trimws(hit[1]), fixed = TRUE)))
  }
  cmd <- Sys.which("sha256sum")
  if (!nzchar(cmd)) stop("sha256sum is required for checksum validation.")
  strsplit(system2(cmd, shQuote(path), stdout = TRUE), "[[:space:]]+")[[1]][1]
}

manifest <- read.csv(manifest_path, stringsAsFactors = FALSE)
missing <- manifest$path[!file.exists(file.path(root, manifest$path))]
if (length(missing)) stop("Missing validation files:\n", paste(missing, collapse = "\n"))
observed <- vapply(file.path(root, manifest$path), sha256_file, character(1))
bad <- which(tolower(manifest$sha256) != observed)
if (length(bad)) {
  stop("Checksum mismatch:\n", paste(manifest$path[bad], collapse = "\n"))
}
cat("Validation passed for ", nrow(manifest), " numerical/data files.\n", sep = "")
