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


moves <- c(1, 0, 1, 0, 0, 0, 1, 2, 1, 0, 1, 0, 2, 2, 0, 1, 1, 0, 0)
start <- c(2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 32, 34, 36, 38)
df <- data.frame(move=move,
                 start=start,
                 non_zero_move=NA,
                 move_2=NA,
                 move_3=NA,
                 start_2=NA,
                 start_3=NA)

# find the non-zero moves
df$non_zero_move <- ifelse(df$move > 0, T, F)
# add an extra start column to add sample index just before a non-zero move
df$start_2 <- ifelse(df$non_zero_move==T, df$start-0.00001, NA)
# add an extra moves column containing original moves shifted by one
df$move_2 <- c(0, df$move[1:NROW(df$move)-1])
# add NAs to useless moves in the extra moves column
df$move_2 <- ifelse(!(is.na(df$start_2)), df$move_2, NA)
df$start_3 <- ifelse(!(is.na(df$start_2)), df$start_2 + stride, NA)
# now interleave the moves and starts
df1 <- dplyr::select(df, start, move)
df2 <- dplyr::select(df, start_2, move_2)
df2 <- dplyr::rename(df2, start=start_2, move=move_2)
df3 <- dplyr::select(df, start_3, move)
df3 <- dplyr::rename(df3, start=start_3)
df <- gdata::interleave(df2, df1, df3)
# select only the complete cases
df3 <- df[complete.cases(df), ]



library(rbokeh)

df <- data.frame(moves <- c(1, 0, 1, 0, 0, 0, 1, 2, 1, 0, 1, 0, 2, 2, 0, 1, 1, 0, 0),
                 start <- c(2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 32, 34, 36, 38))
p1 <- rbokeh::figure(data=df, width=1000, height=200, title='plot_title')
p1 <- rbokeh::ly_lines(p1, start, move, color='#BF1268', width=3, legend = "y")
p1 <- rbokeh::tool_pan(p1, dimensions = "width")
p1 <- rbokeh::tool_wheel_zoom(p1, dimensions = "width")
print(p1)

p2 <- rbokeh::figure(data=df3, width=1000, height=200, title='plot_title')
p2 <- rbokeh::ly_lines(p2, start, move, color='#BF1268', width=3, legend = "y")
p2 <- rbokeh::tool_pan(p2, dimensions = "width")
p2 <- rbokeh::tool_wheel_zoom(p2, dimensions = "width")
print(p2)
p <- rbokeh::grid_plot(list(p1, p2), nrow = 2, link_data = T, same_axes=T)
print(p)



