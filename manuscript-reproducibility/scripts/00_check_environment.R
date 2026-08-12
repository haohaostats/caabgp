root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
required <- c("R", "scripts", "results_r", "case_application")
missing <- required[!file.exists(file.path(root, required))]
if (length(missing)) stop("Run this script from the repository root. Missing: ", paste(missing, collapse = ", "))

if (getRversion() < "4.3.0") warning("R >= 4.3.0 is recommended; found ", getRversion())
cat("R version: ", R.version.string, "\n", sep = "")
cat("Platform: ", R.version$platform, "\n", sep = "")
cat("Logical CPU cores: ", parallel::detectCores(), "\n", sep = "")
cat("pdflatex: ", if (nzchar(Sys.which("pdflatex"))) Sys.which("pdflatex") else "not found (optional)", "\n", sep = "")
cat("Environment check passed.\n")

