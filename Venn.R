#suppressPackageStartupMessages(.libPaths( c( .libPaths(), "/home/clegendre/R/x86_64-unknown-linux-gnu-library/2.15.2") ))

#@@@@@@@@@@@
# FUNCTIONS  
#@@@@@@@@@@@

intersectSeveral <- function(...) { Reduce(intersect, list(...)) } 
getSetListsFromFile<-function(inputFile){
  LF <- read.delim(file=inputFile,header=FALSE)
  return(LF)
  }

#@@@@@@@@@@@
# GET INPUTS
#@@@@@@@@@@@

args<-commandArgs(TRUE)
if (length(args) == 1){
  print("sourcing INI File ...")
  if ( file.exists(args[1] )) { source(args[1])} 
  if (!file.exists(dir_infiles) || dir_infiles == "" ) { stop("dir_iniFiles not DEFINED")} 
}else { stop("INI FILE NOT FOUND") }

#@@@@@@@@@@@
# LOAD LIBS
#@@@@@@@@@@@

suppressPackageStartupMessages(require("VennDiagram"))
suppressPackageStartupMessages(require("gplots"))
suppressPackageStartupMessages(require("Vennerable"))
setwd(dir_infiles) ; cat("Working Directory: ",getwd(),"\n")

#@@@@@@@@@@@
# MAIN
#@@@@@@@@@@@

if (numberOflists == 3){
  ## vennDiagramm 3 lists
  if( !file.exists(file_listA) || !file.exists(file_listB) || !file.exists(file_listC)){stop("ERROR : please provide correct input files (names or path incorrect)")}
  listA<-scan(file = file_listA , what = character(0))
  listB<-scan(file = file_listB , what = character(0))
  listC<-scan(file = file_listC , what = character(0))
  sizeA=length(listA) ; sizeB=length(listB) ; sizeC=length(listC)
  
  print("getting Intersection all lists ...")
  write.table(intersectSeveral(listA,listB,listC), file=paste(prefix_output,".intersections_3L.tsv",sep=""), row.names=F, quote=F, col.names=F)
  print("getting length all lists ...")
  length(intersectSeveral(listA,listB,listC)) ; intersection_l2<-length(intersectSeveral(listA,listB,listC))
  print("getting Union all lists ...")
  union<-(sizeA-intersection_l2)+(sizeB-intersection_l2)+(sizeC-intersection_l2)+intersection_l2
  pct_A=round((sizeA-intersection_l2)/(union)*100, digit=2) ; pct_B=round((sizeB-intersection_l2)/(union)*100, digit=2) ; 
  pct_C=round((sizeC-intersection_l2)/(union)*100, digit=2) ; pct_I=round(intersection_l2/union*100 , digit=2)

  venn.diagram(list(I = listA, II = listB, III = listC),
               category.names = c( category_nameA, category_nameB, category_nameC),
               filename = paste(prefix_output,".Venn.",category_nameA,"_vs_",category_nameB,"_vs_",category_nameC,".png", sep=""),
               lwd = c(3,2,3),
               lty = c(1,1,1),
               fill = c("cornflowerblue", "darkred", "green"), #,"red", "white", "darkblue", "grey",
               label.col = c("black","white","black","white","red","white","black"),
               cex = c(2,2,2,2,4,2,2),
               alpha = 0.65,
               fontface = "bold",
               cat.cex = 0.8,
               cat.fontface = "bold",
               cat.default.pos = "outer",
               cat.pos = c(-45, 45, 180),
               cat.dist = c(0.055, 0.055, 0.055),
               cat.just = list(c(0,0.5),c(0.9,1),c(0.5,-3.5)),
               sub = paste(category_nameA ,"-->", GroupA,"(",sizeA,")  %",pct_A," \t ",category_nameB ,"-->",GroupB," (",sizeB,")  %",pct_B,"  \t", category_nameC ,"-->",GroupC," (",sizeC,")  %",pct_C," \n intersection %",pct_I,"\n",comments2add,"\n"),
               sub.cex = 0.5, sub.col = "darkblue", sub.fontface = "bold",
               main = c(main_title),
               #main.pos = c(0.1,0),
               main.col = "darkblue",
               main.cex = 1.5,
               main.just = c(1, 10),
               sp.cases = TRUE,
               rotation = 1,
               euler.d = F,
               scaled = F,
               )

              if((sizeA>0) && (sizeB>0) && (sizeC>0) ){
               ## using Vennerable to calculate the intersections
               print("Processing Vennerable functions ...")
               VennList<-list(A = listA, B = listB, C = listC)
               names(VennList)<-c(category_nameA,category_nameB,category_nameC)
               V<-VennFromSets(VennList)
               pdf(paste(prefix_output,".Venn.weighted.3L.pdf",sep="")) ; plot(V, doWeights = TRUE, type="circles"); dev.off()
               numOfCombinations<-length(V@IntersectionSets)
               Fname_inter=paste(prefix_output,"sites_in_each_section_sets.3L.tsv",sep=".")
               if(file.exists(Fname_inter)) { file.remove(Fname_inter) }
                  for(i in 2:length(V@IntersectionSets)){
                     if(length(V@IntersectionSets[i][[1]])>0)
                       write.table(t(as.matrix(V@IntersectionSets[i][[1]])), file=Fname_inter, quote=FALSE, row.names=names(V@IntersectionSets[i]), col.names=FALSE, sep="\t", append=TRUE)
               }
              }
}

if (numberOflists == 2){
  if( !file.exists(file_listA) || !file.exists(file_listB)){stop("ERROR : please provide correct input files (names or path incorrect)")}  
  listA<-scan(file = file_listA , what = character(0))
  listB<-scan(file = file_listB , what = character(0))
  sizeA=length(listA) ; sizeB=length(listB)
  
  print("getting Intersection lists and percentages ...")
  write.table(intersectSeveral(listA,listB), file=paste(prefix_output,".intersections_2L.tsv",sep=""), row.names=F, quote=F, col.names=F)
  length(intersectSeveral(listA,listB)) ; intersection_l2<-length(intersectSeveral(listA,listB))
  union<-(sizeA-intersection_l2)+(sizeB-intersection_l2)+intersection_l2
  pct_A=round((sizeA-intersection_l2)/(union)*100, digit=2) ; 
  pct_B=round((sizeB-intersection_l2)/(union)*100, digit=2) ;
  pct_I=round(intersection_l2/union*100 , digit=2)
  
  #venn Diagram for 2 lists
  venn.diagram(list(A = listA, B = listB), 
               paste(prefix_output,".Venn.",category_nameA,"_vs_",category_nameB,".png", sep=""),
               category.names = c( category_nameA, category_nameB),
               lwd = c(2,3),
               lty = 1,
               fill = c("cornflowerblue", "darkred"),
               label.col = c("black", "white", "black"),
               cex = 2,
               alpha = c(0.55, 0.60),
               fontface = "bold",
               cat.fontface = "bold",
               cat.dist = c(0.08, 0.03),
               cat.pos = c(-15, 10),
               sub = paste(category_nameA ,"-->", GroupA,"(",sizeA,")  %",pct_A," \t ",
                           category_nameB ,"-->",GroupB," (",sizeB,")  %",pct_B,"  \t",
                           " \n intersection %",pct_I,"\n",comments2add,"\n"),
               sub.cex = 0.5, sub.col = "darkblue", sub.fontface = "bold",
               main = c(main_title),
               #main.pos = c(0.1,0),
               main.col = "darkblue",
               main.cex = 1.25,
               main.just = c(1, 10),
               rotation.degree = 0,
               sep.dist = 20,
               offset = 0.2,
               euler.d = T,
               scaled = F
  )
  if((sizeA>0) && (sizeB>0)){
    ## using Vennerable to calculate the intersections
    print("Processing Vennerable functions ...")
    VennList<-list(A = listA, B = listB)
    names(VennList)<-c(category_nameA,category_nameB)
    V<-VennFromSets(VennList)
    pdf(paste(prefix_output,".Venn.weighted.2L.pdf",sep="")) ; plot(V, doWeights = TRUE, type="circles"); dev.off()
    numOfCombinations<-length(V@IntersectionSets)
    Fname_inter=paste(prefix_output,"sites_in_each_section_sets.2L.tsv",sep=".")
    if(file.exists(Fname_inter)) { file.remove(Fname_inter) }
    for(i in 2:length(V@IntersectionSets)){
      if(length(V@IntersectionSets[i][[1]])>0)
        write.table(t(as.matrix(V@IntersectionSets[i][[1]])), file=Fname_inter, quote=FALSE, row.names=names(V@IntersectionSets[i]), col.names=FALSE, sep="\t", append=TRUE)
    }
  }
  
}

print("Venn Done.")
