library(VennDiagram)
# A
# T<-venn.diagram(list(A=A,B=B,C=C),filename=NULL
#                 ,lwd=1,lty=2,col=c('red','green')
#                 ,fill=c('red','green')
#                 ,cat.col=c('red','green')
#                 ,reverse=TRUE)
# grid.draw(T)
# data_ins_vcf = 1:30727
# data_ins_gz = 1:20160
ins_vcf = 1:30727
ins_gz = c(1:5956,30728:44932)
# data_del_vcf = 1:14194
# data_del_gz = 1:14402
del_vcf = 
del_gz = 
# Length_A<-length(A)
# Length_B<-length(B)
# Length_C<-length(C)
# Length_AB<-length(intersect(A,B))
# Length_BC<-length(intersect(B,C))
# Length_AC<-length(intersect(A,C))
# Length_ABC<-length(intersect(intersect(A,B),C))
T1<-venn.diagram(list(ins_vcf=ins_vcf,ins_gz=ins_gz),filename=NULL
                ,lwd=1,lty=2
                ,col=c('red','green'),fill=c('red','green')
                ,cat.col=c('red','green')
                ,rotation.degree=90)
T2<-venn.diagram(list(del_vcf=del_vcf,del_gz=del_gz),filename=NULL
                 ,lwd=1,lty=2
                 ,col=c('red','green'),fill=c('red','green')
                 ,cat.col=c('red','green')
                 ,rotation.degree=90)
grid.draw(T)
