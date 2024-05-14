data(bouwer)
bw <- as.data.frame(bouwer)

rc <- 4/2/12        # radius of 2 inches
rw <- 8.25/2/12     # radius of screen
Le <- 10            # screen length
y0 <- 1.44          # initial drawdown
Lw <- 17.92         # height of water above screen bottom
wl <- 6.08          # static wl
H  <- 18.92         # height of water level above base
t  <- as.numeric(bouwer$datetime - bouwer$datetime[1])    # elapsed time
y  <- wl - bouwer$val                                     # change in wl


a <- bouwer_rice(t, y, rw, rc, Le, Lw, H)

expect_equivalent(a * 86400, 4.5, tolerance = 1e-2)


#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
# frec1 = recipe(formula = formula, data = dat) |>
#   step_slug_bouwer_rice() |>
#   plate("dt")
#
# frec2 = Recipe$new(formula = formula, data = dat)$
#   add_step(StepSlugBouwerRice$new())$
#   plate("dt")
#
#
# expect_equivalent(frec1, frec2,
#                             info = "R6 and hydrorecipes api are equivalent")


#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

