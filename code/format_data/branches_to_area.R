library(readxl)
library(dplyr)

d <- read.csv('C:/Lupine length width branches.csv') %>% 
        mutate( area = (LENGTH/2) * (WIDTH/2) * pi ) %>% 
        mutate( NUMBRANCHES2 = NUMBRANCHES^2 )

mod1 <- lm(area ~ NUMBRANCHES, data=d)
mod2 <- lm(area ~ NUMBRANCHES + NUMBRANCHES2, data=d)

xSeq <- seq(min(d$NUMBRANCHES),max(d$NUMBRANCHES), length.out = 100)
ySeq <- coef(mod2)[1] + coef(mod2)[2]*xSeq + coef(mod2)[3]*xSeq^2

#
tiff('C:/brances_to_area_lupine.tiff', unit="in", width=6.3, height=6.3, res=600,compression="lzw")
par(mar=c(4,4, 0.5,0.5))
plot(area ~ NUMBRANCHES, data=d, 
     col = as.factor(POPULATION), pch = 16)
lines(xSeq,ySeq, lwd=2)
legend('topleft', bty='n',legend = unique(d$POPULATION), pch = 16, 
       col = unique(d$POPULATION) )
dev.off()
