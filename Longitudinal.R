##oct23 create longtitudinal data review
darpa=readxl::read_xlsx('/Users/xiaoyi/Documents/2018Fall/Clinical/DARPA Phase 1.xlsx',na = 'NA',sheet = 'All TPs',skip = 4)
colnames(darpa)
colnames(darpa)[7:10]=c('BDNF','IFN-gamma','IL-6','TNF-alpha')


setwd('/Users/xiaoyi/Documents/2018Fall/Clinical/')
dir.create('longitudinal review')
setwd('longitudinal review')

#######function for plot 
lgplot<-function(gr,trt,biom){
  gp=subset(darpa,(Group==gr&Treatment==trt))
  gpname=paste(gr,trt)
pdf(paste0(biom,'_',gpname,'.pdf'))
plot (gp[[biom]] ~gp[["Time Point"]], xlab = "Timepoint", ylab = paste0(biom,"level",sep=" "), type ="n",xaxt="n" )
title (main =gpname, cex = 1.5)
for (i in unique(gp$`Subject ID`)){
  y =gp[[biom]][gp$`Subject ID` == i]; x =gp$`Time Point` [gp$`Subject ID`==i]
  lines (y~x, type = "b", col = i, cex = .75, lwd = .5)
  axis(1,at=as.numeric(unique(gp$`Time Point`)),label=unique(gp$`Time Point`))
}
dev.off()
}

lgplot('PTSD', 'active','BDNF')



for( i in c('BDNF','IFN-gamma','IL-6','TNF-alpha')){
  for (j in c('PTSD','non-PTSD')){
    for (h in c('active','sham')){
      lgplot(j,h,i)
    }
    }
}
