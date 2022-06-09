.onLoad <- function(libname, pkgname) {
  assignInNamespace('tidy.step_ns', value = hydrorecipes:::tidy.step_ns, ns = 'recipes')
  assignInNamespace('tidy.recipe', value = hydrorecipes:::tidy.recipe, ns = 'recipes')
  assignInNamespace('tidy.step', value = hydrorecipes:::tidy.step, ns = 'recipes')
  assignInNamespace('tidy.check', value = hydrorecipes:::tidy.check, ns = 'recipes')
}
