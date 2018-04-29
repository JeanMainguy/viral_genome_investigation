library(ggplot2)
data = read.csv(file = "../result/polyprot_refseq_overview.csv", sep = '\t', header = TRUE, stringsAsFactors = FALSE)

data$polyprotein = ifelse(data$annotation != "None", "Polyprotein" , "No Polyprotein")

data[data == "Genomoviridae"] = "ssDNA viruses"
data[data == "unclassified archaeal viruses" | data == "unclassified bacterial viruses" ] = "unclassified viruses"
data[data == "unclassified virophages" | data == "unclassified viroids" ] = "unclassified viruses"




nodes = table(data$node1)
nodes = data.frame(nodes, stringsAsFactors = FALSE)

colnames(nodes) = c('Node', 'Freq')
other_nodes = nodes$Node[nodes$Freq < 38]
data$node1[data$node1 %in% other_nodes] = "Other"


table = data.frame(table(data$annotation, data$node1))
colnames(table) = c('annotation', 'node1', 'nb')

table = table[order(table$nb), ]



p<-ggplot(data=table, aes(x=reorder(node1, -nb), y=nb, fill=reorder(annotation, nb) )) +geom_bar(stat="identity") + #,  position=position_dodge()) +
  #geom_text(aes(label=Nb_domains), vjust=1.6, color="black",position = position_dodge(0.9), size=3.5)+
  theme_minimal() + theme(axis.text.x = element_text(angle = 90, size=15, vjust=0.5, hjust=1),
                          axis.text.y = element_text(size=12), axis.title=element_text(size=13,face="bold"), 
                          legend.text = element_text( size=12), legend.title = element_text(size=12, face='bold')) + 
 scale_fill_manual(values=c("#44b006","#C70039",  "#999999"))+
  labs( y="Number of RefSeq viral genomes", fill="Polyprotein:", x=NULL)
  #scale_fill_discrete(breaks=c("annotated", "unannotated", 'no_polyprotein'),  labels=c("Annotated", "Not annotated", "Absence" ), color=c("#44b006","#C70039",  "#999999"))

png(filename="../Figures/Overview_RefSeq_Polyprotein_annotation",  width = 700, height = 500)
p
dev.off()
  

  
  
  table = data.frame(table(data$polyprotein, data$node1))
  colnames(table) = c('annotation', 'node1', 'nb')
  
  table = table[order(table$nb), ]
  
p<-ggplot(data=table, aes(x=reorder(node1, -nb), y=nb, fill=reorder(annotation, nb) )) +geom_bar(stat="identity") + #,  position=position_dodge()) +
    #geom_text(aes(label=Nb_domains), vjust=1.6, color="black",position = position_dodge(0.9), size=3.5)+
    theme_minimal() + theme(axis.text.x = element_text(angle = 90, size=15, vjust=0.5, hjust=1),
                            axis.text.y = element_text(size=12), axis.title=element_text(size=13,face="bold"), 
                            legend.text = element_text( size=12), legend.title = element_text(size=12, face='bold')) + 
    scale_fill_manual(values=c("#e5cc3b",  "#999999"))+ 
    labs( y="Number of RefSeq viral genomes", fill=NULL, x=NULL)
  #scale_fill_discrete(breaks=c("annotated", "unannotated", 'no_polyprotein'),  labels=c("Annotated", "Not annotated", "Absence" ), color=c("#44b006","#C70039",  "#999999"))
  
png(filename="../Figures/Overview_RefSeq_Polyproteins",  width = 700, height = 500)
p
dev.off()

