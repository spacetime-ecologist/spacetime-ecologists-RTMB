# Custom function for adding legend in multi-panel figure
# https://stackoverflow.com/questions/52975447/reorganize-sf-multi-plot-and-add-a-legend
add_legend <-
    function(legend,
             col = sf.colors(),
             legend_x = c(0.9, 1.0),
             legend_y = c(0.05, 0.45),
             text_col = "black",
             ...) {
        # Get the axis limits and calculate size
        axisLimits <- par()$usr
        xLength <- axisLimits[2] - axisLimits[1]
        yLength <- axisLimits[4] - axisLimits[3]

        xl <- (1 - legend_x[1]) * par("usr")[1] + (legend_x[1]) * par("usr")[2]
        xr <- (1 - legend_x[2]) * par("usr")[1] + (legend_x[2]) * par("usr")[2]
        yb <- (1 - legend_y[1]) * par("usr")[3] + (legend_y[1]) * par("usr")[4]
        yt <- (1 - legend_y[2]) * par("usr")[3] + (legend_y[2]) * par("usr")[4]
        if (diff(legend_y) > diff(legend_x)) {
            align <- c("lt", "rb")[2]
            gradient <- c("x", "y")[2]
        } else {
            align <- c("lt", "rb")[1]
            gradient <- c("x", "y")[1]
        }

        # Add the legend
        plotrix::color.legend(
            xl = xl,
            xr = xr,
            yb = yb,
            yt = yt,
            legend = legend,
            rect.col = col,
            gradient = "y",
            col = text_col,
            ...
        )
    }
