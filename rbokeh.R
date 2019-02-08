library(rbokeh)

df <- data.frame(x=c(1, 1.999, 2, 2.999, 3, 3.999, 4, 4.999, 5, 5.999),
                 y=c(1, 1, 0, 0, 1, 1, 0, 0, 1, 1))
p1 <- rbokeh::figure(data=df, width=1000, height=200, title='plot_title')
p1 <- rbokeh::ly_lines(p1, x, y, color='#BF1268', width=3, legend = "y")
p1 <- rbokeh::tool_pan(p1, dimensions = "width")
p1 <- rbokeh::tool_wheel_zoom(p1, dimensions = "width")

df2 <- data.frame(x=c(1, 2, 3, 4, 5),
                 y=c(1, 0, 1, 0,  1))
p2 <- rbokeh::figure(data=df2, width=1000, height=200, title='plot_title')
p2 <- rbokeh::ly_lines(p2, x, y, color='#BF1268', width=3, legend = "y")
p2 <- rbokeh::tool_pan(p2, dimensions = "width")
p2 <- rbokeh::tool_wheel_zoom(p2, dimensions = "width")

p <- rbokeh::grid_plot(list(p1, p2), nrow = 2, link_data = T, same_axes=T)
print(p)







