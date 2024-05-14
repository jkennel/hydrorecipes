formula <- as.formula(y~.)

dat <- data.frame(x = rnorm(200),
                  y = rnorm(200))

formula <- as.formula(.~x+y)
frec1 = hydrorecipes:::Recipe$new(formula = formula, data = dat)$
  add_step(hydrorecipes:::StepPgram$new(c(x, y)))$
  plate("df")

frec2 = recipe(formula = formula, data = dat) |>
  step_fft_pgram(c(x,y)) |>
  plate("df")

expect_equivalent(frec1, frec2,
                            info = "R6 and hydrorecipes api are equivalent")


n <- 10000
x <- rnorm(n)
y <- x * 0.2
m <- matrix(c(y,x), ncol = 2)
gain <- as.numeric(Mod(hydrorecipes:::transfer_pgram(m,
                                                 spans = 3,
                                                 detrend = FALSE,
                                                 demean = FALSE,
                                                 taper = 0.1)))

expect_equivalent(gain, rep(0.2, n),
                            info = "transfer_pgram gives the right gain")


n_groups <- 50
gain = as.numeric(Mod(hydrorecipes:::transfer_pgram_smooth(m,
                                                       spans = 3,
                                                       detrend = FALSE,
                                                       demean = FALSE,
                                                       taper = 0.1,
                                                       power = 3,
                                                       n_groups = n_groups)))
expect_equivalent(gain, rep(0.2, n_groups),
                            info = "transfer_pgram_smooth gives the right gain")

