
# Quilty & Roeloffs 1991 --------------------------------------------------

# tf <- fread('/media/kennel/Data/tmp/turkey_flat.csv')
# tf <- fread('/media/kennel/Data/tmp/stockdale_mountain.csv')


# Bouwer & Rice 1976 ------------------------------------------------------
#' @title
#' bouwer
#'
#' @description Keeth elementary dataset 1997-08-14
#' (DATA sheet Slug_Bouwer-Rice.xls)
#' https://pubs.usgs.gov/of/2002/ofr02197/spreadsheets/Slug_Bouwer-Rice.xls
#' @docType data
#'
#' @usage data(bouwer)
#'
#' @format data.table
#'
#' @keywords datasets
#'
#'
#' @examples
#' data(bouwer)
"bouwer"

#' @title
#' bouwer_1976_abc
#'
#' @description a, b and c values for Bouwer and Rice, 1976 solution, for testing
#' https://pubs.usgs.gov/of/2002/ofr02197/spreadsheets/Slug_Bouwer-Rice.xls
#' @docType data
#'
#' @usage data(bouwer_1976_abc)
#'
#' @format data.table
#'
#' @keywords datasets
#'
#'
#' @examples
#' data(bouwer_1976_abc)
"bouwer_1976_abc"


#' @title
#' bouwer_1976_water_level
#'
#' @description dataset for testing
#' https://pubs.usgs.gov/of/2002/ofr02197/spreadsheets/Slug_Bouwer-Rice.xls
#'
#' @docType data
#'
#' @usage data(bouwer_1976_water_level)
#'
#' @format data.table
#' \describe{
#'   \item{datetime}{date and time}
#'   \item{val}{depth to water in feet}
#' }
#' @keywords datasets
#'
#'
#' @examples
#' data(bouwer_1976_water_level)
"bouwer_1976_water_level"



# Rojstaczer & Riley 1990 -------------------------------------------------


#' @title Rojstaczer and Riley (1990) Figure 2 Digitized
#'
#' @description Amplitude of water table response to Earth tides as a function
#' of \eqn{\Omega'}
#'
#' @format A \code{data.table} The columns are:
#' \describe{
#'  \item{\code{ohm}}{water table parameter}
#'  \item{\code{response}}{amplitude of the response}
#' }
#'
#' @examples
#' utils::data(rojstaczer_1990_fig_2)
#' plot(response~ohm, rojstaczer_1990_fig_2, log = 'x')
'rojstaczer_1990_fig_2'



#' @title Rojstaczer and Riley (1990) Figure 3 Digitized
#'
#' @description Areal strain sensitivity and phase of phreatic well to Earth
#' tides as a function of dimensionless frequency Qu and partial penetration.
#' Strain sensitivity is 0.5 (I think this is wrong in the paper where it is
#' listed as 0.05)
#'
#' @format A \code{data.table} The columns are:
#' \describe{
#'  \item{\code{Qu}}{dimensionless frequency}
#'  \item{\code{response}}{gain and phase of response}
#'  \item{\code{b_div_d}}{partial penetration saturated
#'    thickness / aquifer thickness}
#'  \item{\code{variable}}{gain or phase}
#' }
#'
#' @examples
#' utils::data(rojstaczer_1990_fig_3)
'rojstaczer_1990_fig_3'



#' @title Rojstaczer and Riley (1990) Figure 4 Digitized
#'
#' @description Barometric efficiency and phase of phreatic well to atmospheric
#' loading as a function of dimensionless frequency Qu and R/Qu and fully
#' penetrating. Static barometic efficiency is 0.5.
#'
#' @format A \code{data.table} The columns are:
#' \describe{
#'  \item{\code{Qu}}{dimensionless frequency}
#'  \item{\code{response}}{gain and phase of response}
#'  \item{\code{R_div_Qu}}{ratio comparing water table drainage and unsaturated zone influences}
#'  \item{\code{variable}}{gain or phase}
#' }
#'
#' @examples
#' utils::data(rojstaczer_1990_fig_4)
'rojstaczer_1990_fig_4'




# Rojstaczer 1988 Fluid Flow Properties -----------------------------------
#' @title Rojstaczer (1988b) Figure 3 Digitized
#'
#' @description Barometric efficiency and phase as a function of Q/W. Static
#' confined barometric efficiency is 0.5 and R << Q. S and S' are 0.0001.
#'
#' @format A \code{data.table} The columns are:
#' \describe{
#'  \item{\code{W}}{dimensionless frequency}
#'  \item{\code{response}}{gain and phase of response}
#'  \item{\code{Q_div_W}}{ratio comparing water table drainage and borehole storage}
#'  \item{\code{variable}}{gain or phase}
#' }
#'
#' @examples
#' utils::data(rojstaczer_1988b_fig_3)
'rojstaczer_1988b_fig_3'


#' @title Rojstaczer (1988b) Figure 5 Digitized
#'
#' @description Amplitude and phase response as a function of S.
#'
#' @format A \code{data.table} The columns are:
#' \describe{
#'  \item{\code{W}}{dimensionless frequency}
#'  \item{\code{response}}{gain and phase of response}
#'  \item{\code{S}}{storativity}
#'  \item{\code{variable}}{gain or phase}
#' }
#'
#' @examples
#' utils::data(rojstaczer_1988b_fig_5)
'rojstaczer_1988b_fig_5'


# Hsieh 1987 --------------------------------------------------------------

#' @title Hsieh (1987) Figures 2 and 3 Digitized
#'
#' @description Amplitude and phase response as a function of S.
#'
#' @format A \code{data.table} The columns are:
#' \describe{
#'  \item{\code{dimensionless_frequency}}{dimensionless frequency}
#'  \item{\code{response}}{gain and phase of response}
#'  \item{\code{S}}{storativity}
#'  \item{\code{variable}}{gain or phase}
#' }
#'
#' @examples
#' utils::data(hsieh_1987_fig_2_3)
'hsieh_1987_fig_2_3'



# Liu 1989 ----------------------------------------------------------------

#' @title Liu (1989) Figure 8 Digitized
#'
#' @description Amplitude comparison of Cooper and Liu.
#'
#' @format A \code{data.table} The columns are:
#' \describe{
#'  \item{\code{period}}{dimensionless frequency}
#'  \item{\code{response}}{gain of response}
#'  \item{\code{aquifer_thickness}}{thickness of the aquifer (m)}
#'  \item{\code{method}}{cooper or liu method}
#' }
#'
#' @examples
#' utils::data(liu_1989_fig_8)
'liu_1989_fig_8'





# rojstaczer_1990_fig_2 <- fread('/media/kennel/Data/tmp/roj_1990_fig2.csv')
# setnames(rojstaczer_1990_fig_2, c('ohm', 'response'))
# plot(response~ohm, rojstaczer_1990_fig_2, log = 'x')
# use_data(rojstaczer_1990_fig_2)



# g1 <- fread('/media/kennel/Data/tmp/roj1990_f3g_1.csv')
# g1[, b_div_d := 1]
# g2 <- fread('/media/kennel/Data/tmp/roj1990_f3g_50.csv')
# g2[, b_div_d := 0.5]
# g3 <- fread('/media/kennel/Data/tmp/roj1990_f3g_25.csv')
# g3[, b_div_d := 0.25]
# g4 <- fread('/media/kennel/Data/tmp/roj1990_f3g_01.csv')
# g4[, b_div_d := 0.01]
#
# g <- rbindlist(list(g1, g2, g3, g4))
# g[, variable := 'gain']
# setnames(g, c('Qu', 'response', 'b_div_d', 'variable'))
#
# p1 <- fread('/media/kennel/Data/tmp/roj1990_f3p_1.csv')
# p1[, b_div_d := 1]
# p2 <- fread('/media/kennel/Data/tmp/roj1990_f3p_50.csv')
# p2[, b_div_d := 0.5]
# p3 <- fread('/media/kennel/Data/tmp/roj1990_f3p_25.csv')
# p3[, b_div_d := 0.25]
# p4 <- fread('/media/kennel/Data/tmp/roj1990_f3p_01.csv')
# p4[, b_div_d := 0.01]
#
# p <- rbindlist(list(p1, p2, p3, p4))
# p[, variable := 'phase']
# setnames(p, c('Qu', 'response', 'b_div_d', 'variable'))
#
# rojstaczer_1990_fig_3 <- rbind(g, p)
# library(ggplot2)
# ggplot(roj_1990_fig_3, aes(x = Qu, y = response)) + geom_line(aes(color = as.factor(b_div_d))) + scale_x_log10() + facet_wrap(~variable, scales = 'free_y')
# use_data(rojstaczer_1990_fig_3)


#
# g1 <- fread('/media/kennel/Data/tmp/roj1990_f4g_1000.csv')
# g1[, R_div_Qu := 1000]
# g2 <- fread('/media/kennel/Data/tmp/roj1990_f4g_100.csv')
# g2[, R_div_Qu := 100]
# g3 <- fread('/media/kennel/Data/tmp/roj1990_f4g_10.csv')
# g3[, R_div_Qu := 10]
# g4 <- fread('/media/kennel/Data/tmp/roj1990_f4g_1.csv')
# g4[, R_div_Qu := 1]
# g5 <- fread('/media/kennel/Data/tmp/roj1990_f4g_0001.csv')
# g5[, R_div_Qu := 0.0001]
# plot(V2~V1, g1, log = 'x', type='l', ylim = c(0, 1.5))
# points(V2~V1, g2, log = 'x', type='l')
# points(V2~V1, g3, log = 'x', type='l')
# points(V2~V1, g4, log = 'x', type='l')
# points(V2~V1, g5, log = 'x', type='l')
# g <- rbindlist(list(g1, g2, g3, g4, g5))
# g[, variable := 'gain']
# setnames(g, c('Qu', 'response', 'R_div_Qu', 'variable'))
#
# p1 <- fread('/media/kennel/Data/tmp/roj1990_f4p_1000.csv')
# p1[, R_div_Qu := 1000]
# p2 <- fread('/media/kennel/Data/tmp/roj1990_f4p_100.csv')
# p2[, R_div_Qu := 100]
# p3 <- fread('/media/kennel/Data/tmp/roj1990_f4p_10.csv')
# p3[, R_div_Qu := 10]
# p4 <- fread('/media/kennel/Data/tmp/roj1990_f4p_1.csv')
# p4[, R_div_Qu := 1]
# p5 <- fread('/media/kennel/Data/tmp/roj1990_f4p_0001.csv')
# p5[, R_div_Qu := 0.0001]
#
# plot(V2~V1, p1, log = 'x', type='l', ylim = c(-200, -80))
# points(V2~V1, p2, log = 'x', type='l')
# points(V2~V1, p3, log = 'x', type='l')
# points(V2~V1, p4, log = 'x', type='l')
# points(V2~V1, p5, log = 'x', type='l')
#
# p <- rbindlist(list(p1, p2, p3, p4, p5))
# p[, variable := 'phase']
# setnames(p, c('Qu', 'response', 'R_div_Qu', 'variable'))
#
# rojstaczer_1990_fig_4 <- rbind(g, p)
# library(ggplot2)
# ggplot(rojstaczer_1990_fig_4, aes(x = Qu, y = response)) + geom_line(aes(color = as.factor(R_div_Qu))) + scale_x_log10() + facet_wrap(~variable, scales = 'free_y')
# use_data(rojstaczer_1990_fig_4, overwrite = TRUE)



# library(data.table)
# g1 <- fread('/media/kennel/Data/tmp/roj1988b_f3g_10000.csv')
# g1[, Q_div_W := 10000]
# g2 <- fread('/media/kennel/Data/tmp/roj1988b_f3g_1000.csv')
# g2[, Q_div_W := 1000]
# g3 <- fread('/media/kennel/Data/tmp/roj1988b_f3g_100.csv')
# g3[, Q_div_W := 100]
# g4 <- fread('/media/kennel/Data/tmp/roj1988b_f3g_10.csv')
# g4[, Q_div_W := 10]
# g5 <- fread('/media/kennel/Data/tmp/roj1988b_f3g_1.csv')
# g5[, Q_div_W := 1]
#
# plot(V2~V1, g1, log = 'x', type='l', ylim = c(0, 0.8))
# points(V2~V1, g2, log = 'x', type='l')
# points(V2~V1, g3, log = 'x', type='l')
# points(V2~V1, g4, log = 'x', type='l')
# points(V2~V1, g5, log = 'x', type='l')
#
# g <- rbindlist(list(g1, g2, g3, g4, g5))
# g[, variable := 'gain']
# setnames(g, c('W', 'response', 'Q_div_W', 'variable'))
#
# p1 <- fread('/media/kennel/Data/tmp/roj1988b_f3p_10000.csv')
# p1[, Q_div_W := 10000]
# p2 <- fread('/media/kennel/Data/tmp/roj1988b_f3p_1000.csv')
# p2[, Q_div_W := 1000]
# p3 <- fread('/media/kennel/Data/tmp/roj1988b_f3p_100.csv')
# p3[, Q_div_W := 100]
# p4 <- fread('/media/kennel/Data/tmp/roj1988b_f3p_10.csv')
# p4[, Q_div_W := 10]
# p5 <- fread('/media/kennel/Data/tmp/roj1988b_f3p_1.csv')
# p5[, Q_div_W := 1]
#
# plot(V2~V1, p1, log = 'x', type='l', ylim = c(-260, -140))
# points(V2~V1, p2, log = 'x', type='l')
# points(V2~V1, p3, log = 'x', type='l')
# points(V2~V1, p4, log = 'x', type='l')
# points(V2~V1, p5, log = 'x', type='l')
# p <- rbindlist(list(p1, p2, p3, p4, p5))
# p[, variable := 'phase']
# setnames(p, c('W', 'response', 'Q_div_W', 'variable'))
# rojstaczer_1988b_fig_3 <- rbind(g, p)
# library(ggplot2)
# ggplot(rojstaczer_1988b_fig_3, aes(x = W, y = response)) + geom_line(aes(color = as.factor(Q_div_W))) + scale_x_log10() + facet_wrap(~variable, scales = 'free_y')
#
# use_data(rojstaczer_1988b_fig_3, overwrite = TRUE)



# g1 <- fread('/media/kennel/Data/tmp/hsieh_1987_fig3_3.csv')
# g1[, S := 1e-3]
# g2 <- fread('/media/kennel/Data/tmp/hsieh_1987_fig3_4.csv')
# g2[, S := 1e-4]
# g3 <- fread('/media/kennel/Data/tmp/hsieh_1987_fig3_5.csv')
# g3[, S := 1e-5]
# g4 <- fread('/media/kennel/Data/tmp/hsieh_1987_fig3_6.csv')
# g4[, S := 1e-6]
# g5 <- fread('/media/kennel/Data/tmp/hsieh_1987_fig3_7.csv')
# g5[, S := 1e-7]
#
# plot(V2~V1, g1, log = 'x', type='l', ylim = c(0, 1.0))
# points(V2~V1, g2, log = 'x', type='l')
# points(V2~V1, g3, log = 'x', type='l')
# points(V2~V1, g4, log = 'x', type='l')
# points(V2~V1, g5, log = 'x', type='l')
# g <- rbindlist(list(g1, g2, g3, g4, g5))
# g[, variable := 'gain']
# setnames(g, c('dimensionless_frequency', 'response', 'S', 'variable'))
#
#
# p1 <- fread('/media/kennel/Data/tmp/hsieh_1987_fig2_3.csv')
# p1[, S := 1e-3]
# p2 <- fread('/media/kennel/Data/tmp/hsieh_1987_fig2_4.csv')
# p2[, S := 1e-4]
# p3 <- fread('/media/kennel/Data/tmp/hsieh_1987_fig2_5.csv')
# p3[, S := 1e-5]
# p4 <- fread('/media/kennel/Data/tmp/hsieh_1987_fig2_6.csv')
# p4[, S := 1e-6]
# p5 <- fread('/media/kennel/Data/tmp/hsieh_1987_fig2_7.csv')
# p5[, S := 1e-7]
#
# plot(V2~V1, p1, log = 'x', type='l', ylim = c(0, -90))
# points(V2~V1, p2, log = 'x', type='l')
# points(V2~V1, p3, log = 'x', type='l')
# points(V2~V1, p4, log = 'x', type='l')
# points(V2~V1, p5, log = 'x', type='l')
# p <- rbindlist(list(p1, p2, p3, p4, p5))
# p[, variable := 'phase']
# setnames(p, c('dimensionless_frequency', 'response', 'S', 'variable'))
#
# hsieh_1987_fig_2_3 <- rbind(g, p)
# library(ggplot2)
# ggplot(hsieh_1987_fig_2_3, aes(x = dimensionless_frequency, y = response)) + geom_line(aes(color = as.factor(S))) + scale_x_log10() + facet_wrap(~variable, scales = 'free_y')
# use_data(hsieh_1987_fig_2_3, overwrite = TRUE)


#
# library(data.table)
# library(usethis)
# g1 <- fread('/media/kennel/Data/tmp/liu_1989_fig8a_500.csv')
# g1[, aquifer_thickness := 500]
# g2 <- fread('/media/kennel/Data/tmp/liu_1989_fig8a_100.csv')
# g2[, aquifer_thickness := 100]
# g3 <- fread('/media/kennel/Data/tmp/liu_1989_fig8a_50.csv')
# g3[, aquifer_thickness := 50]
# g4 <- fread('/media/kennel/Data/tmp/liu_1989_fig8a_10.csv')
# g4[, aquifer_thickness := 10]
# plot(V2~V1, g1, type='l')
# points(V2~V1, g2, type='l')
# points(V2~V1, g3, type='l')
# points(V2~V1, g4, type='l')
# g <- rbindlist(list(g1, g2, g3, g4))
# setnames(g, c('period', 'response', 'aquifer_thickness'))
# g[, method := 'cooper']
#
#
# p1 <- fread('/media/kennel/Data/tmp/liu_1989_fig8b_500.csv')
# p1[, aquifer_thickness := 500]
# p2 <- fread('/media/kennel/Data/tmp/liu_1989_fig8b_100.csv')
# p2[, aquifer_thickness := 100]
# p3 <- fread('/media/kennel/Data/tmp/liu_1989_fig8b_50.csv')
# p3[, aquifer_thickness := 50]
# p4 <- fread('/media/kennel/Data/tmp/liu_1989_fig8b_10.csv')
# p4[, aquifer_thickness := 10]
#
# plot(V2~V1, p1, type='l')
# plot(V2~V1, p2, type='l')
# plot(V2~V1, p3, type='l')
# plot(V2~V1, p4, type='l')
#
# p <- rbindlist(list(p1, p2, p3, p4))
# setnames(p, c('period', 'response', 'aquifer_thickness'))
# p[, method := 'liu']
# liu_1989_fig_8 <- rbind(p, g)
# setkey(liu_1989_fig_8, method, aquifer_thickness, period)
#
# use_data(liu_1989_fig_8, overwrite = TRUE)
#
# ggplot(liu_1989_fig_8, aes(x = period, y = response)) + geom_line(aes(color = as.factor(aquifer_thickness))) + scale_x_log10() + facet_wrap(~method, scales = 'free_y')
#
# plot(response~period, liu_1989_fig_8[method =='cooper' & aquifer_thickness == 10], type='o', xlim = c(20, 21))






# g1 <- fread('/media/kennel/Data/tmp/roj1988b_f5g_7.csv')
# g1[, S := 1e-7]
# g2 <- fread('/media/kennel/Data/tmp/roj1988b_f5g_5.csv')
# g2[, S := 1e-5]
# g3 <- fread('/media/kennel/Data/tmp/roj1988b_f5g_3.csv')
# g3[, S := 1e-3]
#
# plot(V2~V1, g1, type='l', log = 'x')
# points(V2~V1, g2, type='l')
# points(V2~V1, g3, type='l')
#
# g <- rbindlist(list(g1, g2, g3))
# setnames(g, c('W', 'response', 'S'))
# g[, variable := 'gain']
#
#
#
# p1 <- fread('/media/kennel/Data/tmp/roj1988b_f5p_7.csv')
# p1[, S := 1e-7]
# p2 <- fread('/media/kennel/Data/tmp/roj1988b_f5p_5.csv')
# p2[, S := 1e-5]
# p3 <- fread('/media/kennel/Data/tmp/roj1988b_f5p_3.csv')
# p3[, S := 1e-3]
#
#
# p <- rbindlist(list(p1, p2, p3))
# setnames(p, c('W', 'response', 'S'))
# p[, variable := 'phase']
#
#
#
# rojstaczer_1988b_fig_5 <- rbind(p, g)
#
#
# ggplot(rojstaczer_1988b_fig_5, aes(x = W, y = response)) + geom_line(aes(color = as.factor(S))) + scale_x_log10() + facet_wrap(~variable, scales = 'free_y')
#
# use_data(rojstaczer_1988b_fig_5, overwrite = TRUE)
#
#
#
# plot(V2~V1,   p1, type='l', log = 'x')
# points(V2~V1, p2, type='l')
# points(V2~V1, p3, type='l')




# g1 <- fread('/media/kennel/Data/tmp/roj1988b_f6g_1000.csv')
# g1[, R_div_Q := 1000]
# g2 <- fread('/media/kennel/Data/tmp/roj1988b_f6g_100.csv')
# g2[, R_div_Q := 100]
# g3 <- fread('/media/kennel/Data/tmp/roj1988b_f6g_10.csv')
# g3[, R_div_Q := 10]
# g4 <- fread('/media/kennel/Data/tmp/roj1988b_f6g_1.csv')
# g4[, R_div_Q := 1]
# g5 <- fread('/media/kennel/Data/tmp/roj1988b_f6g_0001.csv')
# g5[, R_div_Q := 0.0001]
#
#
# plot(V2~V1, g1, type='l', log = 'x', ylim = c(0, 1.5))
# points(V2~V1, g2, type='l')
# points(V2~V1, g3, type='l')
# points(V2~V1, g4, type='l')
# points(V2~V1, g5, type='l')
#
# g <- rbindlist(list(g1, g2, g3, g4, g5))
# g[, variable := 'gain']
#
# setnames(g, c('dimensionless_frequency', 'response', 'R_div_Q', 'variable'))
#
#
# p1 <- fread('/media/kennel/Data/tmp/roj1988b_f6p_1000.csv')
# p1[, R_div_Q := 1000]
# p2 <- fread('/media/kennel/Data/tmp/roj1988b_f6p_100.csv')
# p2[, R_div_Q := 100]
# p3 <- fread('/media/kennel/Data/tmp/roj1988b_f6p_10.csv')
# p3[, R_div_Q := 10]
# p4 <- fread('/media/kennel/Data/tmp/roj1988b_f6p_1.csv')
# p4[, R_div_Q := 1]
# p5 <- fread('/media/kennel/Data/tmp/roj1988b_f6p_0001.csv')
# p5[, R_div_Q := 0.0001]
#
#
# plot(V2~V1, p1, type='l', log = 'x', ylim = c(-200, -100))
# points(V2~V1, p2, type='l')
# points(V2~V1, p3, type='l')
# points(V2~V1, p4, type='l')
# points(V2~V1, p5, type='l')
#
# p <- rbindlist(list(p1, p2, p3, p4, p5))
# p[, variable := 'phase']
# setnames(p, c('dimensionless_frequency', 'response', 'R_div_Q', 'variable'))
#
# rojstaczer_1988b_fig_6 <- rbind(p, g)
#
#
# ggplot(rojstaczer_1988b_fig_6, aes(x = dimensionless_frequency, y = response)) + geom_line(aes(color = as.factor(R_div_Q))) + scale_x_log10() + facet_wrap(~variable, scales = 'free_y')
#
# use_data(rojstaczer_1988b_fig_6, overwrite = TRUE)


# g1 <- fread('/media/kennel/Data/tmp/rojstaczer_tidal_gain.csv')
# g1[, variable := 'gain']
# p1 <- fread('/media/kennel/Data/tmp/rojstaczer_tidal_phase.csv')
# p1[, variable := 'phase']
#
# rojstaczer_1988a_fig_3 <- rbind(g1, p1)
# setnames(rojstaczer_1988a_fig_3, c('dimensionless_frequency', 'response', 'variable'))
# use_data(rojstaczer_1988a_fig_3, overwrite = TRUE)
