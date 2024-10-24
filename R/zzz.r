.onAttach <- function(...) {
  info <- c("New version of synthpop with disclosure functions",
            "Find out more at https://www.synthpop.org.uk/")
    packageStartupMessage(paste(strwrap(info), collapse = "\n"))
}
