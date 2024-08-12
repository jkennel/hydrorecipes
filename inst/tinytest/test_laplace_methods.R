n <- 1000000
time = c(1:n)
Tr = 200.0 # transmissivity of aquifer, m^2/d
S = 0.0005 # storage coefficient of aquifer, -
Q = 800.0  # discharge of well, m^3/d
rw = 10.0  # radius to well, m
n_terms <- 18L
prec = 1e-5

a <- hydrorecipes:::theis_laplace(time, rw, Tr, S, Q, prec, n_terms)
b <- hydrorecipes::theis_aniso_time(0, rw, S, Tr, Tr, 1, time, rep(Q, n))[[1]]

bench::mark(a <- hydrorecipes:::theis_laplace(time, rw, Tr, S, Q, prec, n_terms),
            b <- hydrorecipes::theis_aniso_time(0, rw, S, Tr, Tr, 1, time, rep(Q, n))[[1]],
            check = FALSE)

expect_equivalent(a, b, info = "theis laplace and standard method are equal", tolerance = 1e-6)


cc <- 5e4
lab = sqrt(cc * Tr)
rw = 2500 # radius to well, m
lab = sqrt(cc * Tr)

a <- hydrorecipes:::hantush_jacob_laplace(time, cc, rw, Tr, S, Q, prec, n_terms)
b <- hydrorecipes::hantush_jacob(time, rep(Q, n), rw, S, Tr, lab, prec)[[1]]

# bench::mark(a <- hydrorecipes:::hantush_jacob_laplace(time, cc, rw, Tr, S, Q, prec, n_terms),
#             b <- hydrorecipes::hantush_jacob(time, rep(Q, n), rw, S, Tr, lab, prec)[[1]],
#             check = FALSE)

expect_equivalent(a, b, info = "hantush-jacob laplace and standard method are equal", tolerance = 1e-6)

