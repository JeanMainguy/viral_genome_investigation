library(ggplot2)
library(scales)

data = read.csv("../result/polyproteins_interpro_domains.csv", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

overlap_data = data[data$overflow == 'True', ]
overlap_data = overlap_data[overlap_data$Expected_peptides_number == 'True', ]


#Group overlap distance
overlap_data$groupedby5 <-cut(overlap_data$overflow_distance , seq(0, 150, 10))

#DISTANCE OF OVERLAPING
overlap_database = data.frame( table(overlap_data$groupedby5, overlap_data$app))
colnames(overlap_database) = c("overlap", "method", "Nb_domains")

p<-ggplot(data=overlap_database, aes(x=overlap, y=Nb_domains)) +
  geom_bar(stat="identity",fill = "#f28759") + theme_minimal() + 
  theme(axis.text.x = element_text(angle = 90, size=15, vjust=0.5, hjust=1),
        axis.text.y = element_text(size=15), axis.title=element_text(size=15,face="bold")) +
  labs( y="Number of domain annotations", x="Overlap distance (amino acid)")+
  scale_fill_manual(values=c("red"))

png(filename="../Figures/Overlap_distance_of_domain_annotation",  width = 700, height = 500)
p
dev.off()

#COLOR BY DATABASE
p<-ggplot(data=overlap_database, aes(x=overlap, y=Nb_domains, fill=reorder(method, Nb_domains))) +
  geom_bar(stat="identity") +  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 90, size=15, vjust=0.5, hjust=1),
        axis.text.y = element_text(size=15), axis.title=element_text(size=15,face="bold"), 
        legend.text = element_text( size=15), legend.title = element_text(size=15, face='bold')) +
  labs( y="Number of domain annotations", fill="Interpro Database:", x="Overlap distance (amino acid)")

png(filename="../Figures/Overlap_distance_by_database",  width = 700, height = 500)
p
dev.off()

#COLOR BY DOMAIN ANNOTATIONS
#We don't display the name of domain that appears less than 5 time
#They are label as Other
domains = table(overlap_data$Name)
domains = data.frame(domains, stringsAsFactors = FALSE)
colnames(domains) = c('Domain', 'Freq')
other_domain = domains$Domain[domains$Freq < 5]
overlap_data$Name[overlap_data$Name %in% other_domain] = "Other"

domain_annotations_data = data.frame( table(overlap_data$groupedby5, overlap_data$Name), stringsAsFactors = FALSE)
colnames(domain_annotations_data) = c("overlap", "Domain", "Nb_domains")
domain_annotations_data = rev(domain_annotations_data)

p<-ggplot(data=domain_annotations_data, aes(x=overlap, y=Nb_domains, fill=Domain)) +
  geom_bar(stat="identity") +  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 90, size=15, vjust=0.5, hjust=1),
        axis.text.y = element_text(size=15), axis.title=element_text(size=15,face="bold"), 
        legend.text = element_text( size=15), legend.title = element_text(size=15, face='bold')) +
  labs( y="Number of domain annotations", fill="Domain Name:", x="Overlap distance (amino acid)") +
 scale_fill_manual(values=c("#999999", hue_pal()(21)))
p

png(filename="../Figures/Overlap_distance_by_domain_name",  width = 700, height = 500)
p
dev.off()
