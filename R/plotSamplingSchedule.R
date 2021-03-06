#' Plot sampling schedule
#'
#' This function plots the sampling schedule.
#'
#' @param ps (Required). A \code{\link{phyloseq-class}}.
#' @param time_var Character string. The name of the variable to map time in the plot (integer vector).
#' @param subjectID_var Character string. The name of the variable to map Subject Ids in the plot (factor vector).
#' @param main_factor Character string. The name of the variable to map main factor in the plot (factor vector).
#' @param theme_manual Character string. Specifying the object of \link[ggplot2]{theme_update}.
#'
#' @return \code{\link{ggplot2-class}} that displays the sampling schedule.
#' @export
plotSamplingSchedule <- function(ps,
    time_var,
    subjectID_var,
    main_factor,
    theme_manual = set.theme) {

    samdf <- sample_data(ps) %>% data.frame
    if (!is.integer(samdf[, time_var])) {
        stop(paste(time_var, "must be integer", sep = " "))
    }
    if (!is.character(samdf[, subjectID_var])) {
        stop(paste(subjectID_var, "must be character", sep = " "))
    }
    if (!is.factor(samdf[, main_factor])) {
        stop(paste(main_factor, "must be factor", sep = " "))
    }

    time_brks <- seq(1, max(samdf[, time_var]))

    p <- ggplot(samdf) +
        geom_point(aes_string(x = time_var, y = subjectID_var,
        color = main_factor), position = position_jitter(width = 0.01, height = 0.1), size = 1) +
        theme_set(theme_manual) +
        facet_wrap(as.formula(paste("~",
        main_factor))) +
        scale_x_continuous(breaks = time_brks)

    return(p)
}
