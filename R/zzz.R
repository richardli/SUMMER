.onAttach <- function( lib, pkg ) {
   packageStartupMessage(
      paste0("SUMMER version ", utils::packageDescription("SUMMER")$Version, "\n", 
      		 "  See latest changes with 'news(package = 'SUMMER')'"),
      domain = NULL,  appendLF = TRUE )
}