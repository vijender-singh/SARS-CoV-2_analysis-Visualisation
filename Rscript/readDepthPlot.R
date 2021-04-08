# test if there is at least one argument: if not, return an error

library(ggplot2)
library(tidyr)

if (length(args)==0) {
  stop("Arguments must be supplied (path to readDepth file directory, ouputdirectory).n", call.=FALSE)
} else if (length(args)==1) {
  directory=args[1]
  outputDir=args[2]
}


files=list.files(path=directory, pattern="*Depth.txt")

for (id in seq(1,length(files))){
    if (id ==1){  
      table=read.table(files[id])[,c(2,3)]
      colnames(table)<-c("Cord",strsplit(files[id],"___")[[1]][1])
      }else{
      table1=read.table(files[id])[,c(2,3)]
      colnames(table1)<-c("Cord",strsplit(files[id],"___")[[1]][1])
#      head(table)
#      head(table1)
      table=merge(table,table1,by.x="Cord",by.y="Cord")   
      }
}

table$loc<-table$Cord
table$loc<-factor(table$loc)
data_long <- gather(table, sample, depth, files[1]:files[length(files)], factor_key=TRUE)
data_long[data_long==0] <- 1

ggplot(data_long , aes(x = loc, y = log2(depth))) + 
  geom_line()+ facet_wrap(~sample,  ncol=2, scales = "free")

ggsave(paste0(outputDir,"read_depth.png"))
dev.off()
