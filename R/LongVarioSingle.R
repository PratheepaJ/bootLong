#' Plot the variogram given the taxon index.
#'
#' @param point Logical. Whether variogram with observed values.
#' @inheritParams longCorreloSingle
#' @inheritParams plotSamplingSchedule
#'
#' @return \code{ggplot} object of variogram
#' @export
longVarioSingle <- function(ps, main_factor, time_var, subjectID_var, taxon,
    point = FALSE, taxlevel = "Species") {

    taxon_name = tax_table(ps)[taxon, taxlevel]
    df.taxa = data.frame(sample_data(ps), otu = as.numeric(t(otu_table(ps)[taxon,
        ])))
    names(df.taxa)[names(df.taxa) == main_factor] = "Group"
    names(df.taxa)[names(df.taxa) == time_var] = "Time"
    names(df.taxa)[names(df.taxa) == subjectID_var] = "SubjectID"

    if (!is.numeric(df.taxa$Time)) {
        df.taxa$Time = as.numeric(df.taxa$Time)
    }

    time_brks <- seq(1, max(df.taxa$Time))

    df.taxa.sep = split(df.taxa, df.taxa$Group)

    res.sep = lapply(df.taxa.sep, function(m) {
        joineR::variogram(m$SubjectID, time = m$Time, Y = m$otu)
    })

    df.sep = lapply(res.sep, function(m) {
        data.frame(m$svar)
    })
    Group = names(df.sep)
    sig2 = lapply(res.sep, function(m) {
        m$sigma2
    })
    for (i in 1:length(df.sep)) {
        df.sep[[i]]$Group = Group[i]
        df.sep[[i]]$sigma2 = sig2[[i]]
    }

    sig2 = unlist(sig2)
    df.sep = do.call("rbind", df.sep)
    df.sep$Group = as.factor(df.sep$Group)
    yl <- max(sig2)
    if (point == TRUE) {
        p = ggplot(df.sep) + geom_point(aes(x = vt, y = vv, col = Group, group = Group),
            size = 1) + geom_smooth(aes(x = vt, y = vv, col = Group, group = Group),
            method = "loess") + geom_hline(aes(yintercept = sigma2, group = Group,
            col = Group), linetype = 2, size = 1) + xlab("u") + ylab("v(u)") +
            ggtitle(taxon_name) + theme(plot.title = element_text(hjust = 0.5,
            size = 8), axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8),
            axis.title.x = element_text(size = 8), axis.title.y = element_text(size = 8),
            legend.text = element_text(size = 8), legend.title = element_text(size = 8)) +
            scale_y_continuous(breaks = time_brks)
    } else {
        p = ggplot(df.sep) + geom_point(aes(x = vt, y = vv, col = Group, group = Group),
            size = 1, col = "white") + geom_smooth(aes(x = vt, y = vv, col = Group,
            group = Group), method = "loess") + geom_hline(aes(yintercept = sigma2,
            group = Group, col = Group), linetype = 2, size = 1) + xlab("u") +
            ylab("v(u)") + ggtitle(taxon_name) + theme(plot.title = element_text(hjust = 0.5,
            size = 8), axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8),
            axis.title.x = element_text(size = 8), axis.title.y = element_text(size = 8),
            legend.text = element_text(size = 8), legend.title = element_text(size = 8)) +
            scale_x_continuous(breaks = time_brks)

    }

    return(p)
}
