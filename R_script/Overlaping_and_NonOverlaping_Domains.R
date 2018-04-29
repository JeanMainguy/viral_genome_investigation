library(ggplot2)

data = read.csv("result/polyproteins_interpro_domains.csv", header = TRUE, sep = "\t")
data = data[data$Expected_peptides_number == 'True', ]

#SIMPLE HISTOGRAMME WITH OVERLAPPING AND NON OVERLAPPING
table = data.frame( table(data$included))
colnames(table) = c("Fully_included", "Nb_domains")
table$Overlap = ifelse(table$Fully_included == 'True', 'Non overlapping', 'Overlapping')

p<-ggplot(data=table, aes(x=Overlap, y=Nb_domains, fill=Overlap)) +geom_bar(stat="identity") + theme_minimal() +
  labs(x =NULL, y="Number of domain annotations", fill=NULL) + theme(axis.text=element_text(size=20), axis.title=element_text(size=18,face="bold"), legend.text = element_text( size=18) )

p = p+scale_fill_manual(values=c("#E69F00", "#56B4E9"), guide=FALSE)
p = p+  geom_text(aes(label=Nb_domains), vjust=-0.1, color="black",
                  position = position_dodge(0.95), size=7)

png(filename="Figures/Overlaping_and_NonOverlaping_Domains",  width = 400, height = 500)
p
dev.off()



#Histogram : Overlaping_and_NonOverlaping_Domains_by_database
table = data.frame( table(data$included, data$app))
colnames(table) = c("Fully_included", "Method", "Nb_domains")
table$Method <- factor(table$Method,levels = c( "ProDom", "SMART","CDD","ProSiteProfiles", "Pfam"))

p<-ggplot(data=table, aes(x=Method, y=Nb_domains, fill=Fully_included)) +geom_bar(stat="identity",  position=position_dodge()) +
  theme_minimal() +geom_text(aes(label=Nb_domains), vjust=-0.1, color="black",
                             position = position_dodge(0.95), size=7)+ theme(axis.text=element_text(size=16),
                                                                             axis.title=element_text(size=18,face="bold"), legend.text = element_text( size=18) ) +
  theme(legend.title = element_text(colour="blue", size=20)) +
  scale_fill_discrete(breaks=c("False", "True"),
                      labels=c("Overlaping", "Non Overlaping"))+
  labs(x = "Interproscan Method", y="Number of domain annotations", fill=NULL)

png(filename="Figures/Overlaping_and_NonOverlaping_Domains_by_database",  width = 700, height = 500)
p
dev.off()

# Overlaping_and_NonOverlaping_Domains_by_database_frequence
table_frq = table(data$included, data$app)
table_frq = prop.table(table_frq, 2)

data_frq = data.frame(table_frq)
colnames(data_frq) = c("Fully_included", "Method", "Nb_domains")

data_frq$Nb_domains = round(data_frq$Nb_domains*100, 1)
data_frq$Method <- factor(data_frq$Method,levels = c( "ProDom", "SMART","CDD","ProSite", "Pfam"))

p<-ggplot(data=data_frq, aes(x=Method, y=Nb_domains, fill=Fully_included)) + theme_minimal() +geom_bar(stat="identity",  position=position_dodge()) +
  geom_text(aes(label=Nb_domains), vjust=-0.1, color="black",
            position = position_dodge(0.95), size=7)+ theme(axis.text=element_text(size=20),
                                                            axis.title=element_text(size=18,face="bold"), legend.text = element_text( size=18) ) +
  theme(legend.title = element_text(colour="blue", size=20)) +
  scale_fill_discrete(breaks=c("False", "True"),
                      labels=c("Overlaping", "Non Overlaping"))+

  labs(x = "Interproscan Method", y="Number of domain annotations", fill=NULL)

png(filename="Figures/Overlaping_and_NonOverlaping_Domains_by_database_frequence",  width = 800, height = 500)
p
dev.off()
