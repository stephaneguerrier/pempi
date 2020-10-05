#' @export
gui = function(){
  appDir = system.file("gui", package = "cape")
  shiny::runApp(appDir, display.mode = "normal")
}
