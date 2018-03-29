library(ggplot2)
args = commandArgs(T)

print(args)

input_file = args[1]
output_file = args[2]
ylabel = args[3]
ymax = as.numeric(args[4])
title = args[5]

width=5
height=5

output_png = gsub("pdf", "png", output_file)

D = read.table(input_file, header=T)

pdf(output_file, width=width, height=height)
p = ggplot(data=D, aes(x=SET, y=COUNT, group=TYPE, fill=TYPE)) + 
geom_bar(colour="white", stat="identity",position=position_dodge(),size=0, width=.7) +
scale_fill_manual(name='Set', values=c('ASN'="blue", "SNV"="forestgreen", "UNION"="darkorange2")) + 
coord_cartesian(ylim = c(0, ymax)) + 
ylab(ylabel) + 
xlab("") + 
theme_minimal() + 
theme(axis.text=element_text(size=rel(2)), axis.title=element_text(size=rel(2))) + 
ggtitle(title) 

plot(p)
dev.off()

png(output_png, width=width,height=height,units="in",res=300)
plot(p)
dev.off()
