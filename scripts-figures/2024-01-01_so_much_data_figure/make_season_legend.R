
plot(1:10,1:10, type = "n", axes = F, ann = F)

x.locs <- seq(.25,9.65,along.with = 1:7)
# manual adjust:    seas   ice    spr   clr    earl   late  fall
x.locs <- x.locs + c(.2,   .1,   -.1,   -.4,   -.15,   -.05,   0)

rect.width <- .2
rect.left <- x.locs - rect.width / 2
rect.right <- x.locs + rect.width / 2
text.left <- rect.right + .1
rect.left <- rect.left[-1]
rect.right <- rect.right[-1]

y.locs <- 5
rect.height <- 4.5
rect.top <- y.locs + rect.height / 2
rect.bottom <- y.locs - rect.height / 2

rect(xleft = rect.left, xright = rect.right, ybottom = rect.bottom, ytop = rect.top, 
     col = c(adjustcolor("snow3",.3),
             adjustcolor("tan4",.6),
             adjustcolor("cornflowerblue",.6),
             adjustcolor("chartreuse4",.6),
             adjustcolor("purple",.4),
             adjustcolor("hotpink2",.5)), 
     border = c(adjustcolor("snow3",1),
                adjustcolor("tan4",1),
                adjustcolor("cornflowerblue",1),
                adjustcolor("chartreuse4",1),
                adjustcolor("purple",1),
                adjustcolor("hotpink2",1)))
text(x = text.left, y = y.locs, labels = c("Season:","Ice-On","Spring","Clearwater","Early\nSummer","Late \nSummer","Fall"), adj = 0, xpd = NA)
