data("kennel_2020")
kennel_2020[1e4, wl := 13.36]
frec1 = recipe(wl~baro, data = kennel_2020)$
  add_step(frecipes:::StepCompareColumns$new(data = wl, compare = baro, n_sd = 15))$
  prep()$
  bake()

# sum(frec1$result[[3]])
# which(frec1$result[[3]])
#
# plot(frec1$result[[2]], col = (1 + as.numeric(frec1$result[[3]])),
#      type = "p", pch= 20, cex = 0.2)
