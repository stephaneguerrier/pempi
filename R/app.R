#' @title GUI
#' @description GUI of the cape R package
#' @export
#' @author Mucyo Karemera, Stephane Guerrier
#' @importFrom shiny runApp
gui = function(){
  appDir = system.file("gui", package = "cape")
  runApp(appDir, display.mode = "normal")
}
