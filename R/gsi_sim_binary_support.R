

#' return the path where gsi_sim should be in the R system paths
#' 
#' @keywords internal
gsi_sim_binary_path <- function() {
  file.path(system.file(package = "assigner"), "bin", "gsi_sim")
}

#' return TRUE if gsi_sim exists where it should be
#' @keywords internal
gsi_sim_exists <- function() {
  file.exists(gsi_sim_binary_path())
}


#' return TRUE if gsi_sim is executable
#' @keywords internal
gsi_sim_is_executable <- function() {
  NULL #incomplete
}


#' file path to be used in a call to gsi_sim.
#' 
#' This version checks to make sure it is there and throws an
#' error with a suggestion of how to get it if it is not there.
#' @export
#' @keywords internal
gsi_sim_binary <- function() {
  if (!gsi_sim_exists()) stop("Can't find the gsi_sim executable where it was expected
at ", gsi_sim_binary_path(), ".  
If you have internet access, you can install it
from within R by invoking the function \"install_gsi_sim()\"")
  
  # then I should check to make sure it is executable
  
  # if so, return the path
  gsi_sim_binary_path()
  
}


#' downloads gsi_sim that is appropriate for the operating system
#' 
#' If the system is Mac, or Windows, this function will
#' download a precompiled binary from GitHub.  In other cases, or
#' if fromSource == TRUE, this function will attempt to download
#' the source code and compile the program from source and install 
#' it.
#' 
#' If this function fails, then you can just compile gsi_sim by
#' going to GITHUB_URL and compiling it yourself and naming the
#' executable gsi_sim and putting it at the location specified by the
#' function \code{\link{gsi_sim_binary_path}}.
#' @param commit  The full SHA-1 hash from GitHub from which to get
#' the binary or source
#' @param fromSource If TRUE, download source, even if a binary is available.
#' If FALSE, then it will download a precompiled binary, if available.  If a 
#' binary is not available, then it will attempt to download the source.  
#' @export
# @keywords internal
#' @importFrom utils download.file
#' 
install_gsi_sim <- function(commit = "080f462c8eff035fa3e9f2fdce26c3ac013e208a", fromSource = FALSE) {
  
  # make a bin directory
  suppressWarnings(dir.create(file.path(system.file(package = "assigner"), "bin")))
  
  uname <- Sys.info()["sysname"]
  urlbase <- paste("https://github.com/eriqande/gsi_sim/blob/", commit,
               "/gsi_sim-", sep = "")
  
  if (fromSource == FALSE) {
    if (uname == "Darwin") {
      url <- paste(urlbase, "Darwin", sep = "")
    }
    if (uname == "Windows") {
      url <- paste(urlbase, "MINGW32_NT-6.1", sep = "")
    }
    if (uname == "Darwin" || uname == "Windows") {
      message("Downloading file ", url)
      message("And copying to ", gsi_sim_binary_path())
      utils::download.file(url = url, destfile = gsi_sim_binary_path())
      Sys.chmod(gsi_sim_binary_path()) # make it executable
    }
    return(NULL)
  }
  
  if (uname == "Linux" || fromSource == TRUE) {  # in this case we will just compile from source
    td <- tempdir()
    
    message("Will be cloning gsi_sim repository to ", td)
    message("")
    message("Removing any earlier instances of the repository in that temp directory")
    message("")
    system(paste("cd", td, "; rm -r -f gsi_sim"))
    message("Cloning repository, dealing with submodules, compiling gsi_sim ")
    message("")
    comm <- paste("cd ", td, 
                  "&&  git clone https://github.com/eriqande/gsi_sim.git ",
                  "&&  cd gsi_sim ",
                  "&&  git checkout ", commit,
                  "&& git submodule init && git submodule update ",
                  "&& ./Compile_gsi_sim.sh ", sep = "")
    boing <- system(comm)
    if (boing != 0) {
      stop("Failed trying to clone and compile gsi_sim")
    } else {
      message("Apparently successful compiling gsi_sim.  Now copying to ", gsi_sim_binary_path())
      trycopy <- file.copy(from = paste(td, "/gsi_sim/gsi_sim-", uname, sep = ""), 
                to = gsi_sim_binary_path(),
                overwrite = TRUE)
      if (trycopy == FALSE) stop("Apparently failed trying to copy ", 
                                paste(td, "/gsi_sim/gsi_sim-", uname, sep = ""),
                                "to ",
                                gsi_sim_binary_path())
    }
  }
}



