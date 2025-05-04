#' @include generics.R
NULL
#' @rdname sn_commands
#' @export
setMethod("sn_commands", "Shennong", function(object, ...) {
  return(names(slot(object = object, name = "commands")))
})

