library(scales)


myColors_diagnosis <- c("lightgrey","darkgoldenrod1","dodgerblue3","lightskyblue","lightskyblue","khaki1","khaki1","khaki1","#636b4d","#993399","darkgrey","#819377")
color_possibilities_diagnosis<-c( "Control","UC","CD","IBD-U (CD-like)","ND - CD","IBD-U (UC-like)","UC-like","ND - UC", "IBD-U", "Other.GI","Neonatal","IBD")
names(myColors_diagnosis) <- color_possibilities_diagnosis
fillscale_diagnosis <- scale_fill_manual(name="Diagnosis",
                               values = myColors_diagnosis, drop = T, limits=force)
colscale_diagnosis <- scale_color_manual(name="Diagnosis",
                                         values = myColors_diagnosis, drop = T, limits=force)





myColors_sampsite <- c("#1a9850","#a6d96a","cornflowerblue","cornflowerblue","#a6d96a","grey","#a6d96a","cornflowerblue","#F75C03","#F75C03")
color_possibilities_sampsite<-c( "AC","SC","TI","proximal","distal","other","sigmoid colon","terminal ileum","Duo","DUO")
names(myColors_sampsite) <- color_possibilities_sampsite
fillscale_sampsite <- scale_fill_manual(name="Sample Site",
                                         values = myColors_sampsite, drop = T, limits=force)
colscale_sampsite <- scale_color_manual(name="Sample Site",
                                        values = myColors_sampsite, drop = T, limits=force)





th <-   theme(axis.text=element_text(size=10),
              axis.title=element_text(size=12),
              strip.text = element_text(size = 12),
              legend.text=element_text(size=12),
              legend.title=element_text(size=14))

th_present <- theme(axis.text=element_text(size=12),
                    axis.title=element_text(size=14),
                    strip.text = element_text(size = 12),
                    legend.text=element_text(size=12),
                    legend.title=element_text(size=14))



pass_col<-c("grey30","#9E0142", "#D53E4F", "#F46D43", "#FDAE61","#FEC776","#FEE08B", "#FFFFBF","#E6F598", "#C9E99E","#ABDDA4","#66C2A5","#4CA5B1","#3288BD", "#5E4FA2","#762a83","#3f007d")
names(pass_col)<-c(0,1,2,3,4,5,6,7,8,9,10,11,12,14,16,21,23)



myColors_diagnosis_sample_samplesize <- c("lightgrey","grey80","darkgoldenrod1","darkgoldenrod1","dodgerblue3","#74abe1","#819377")
names(myColors_diagnosis_sample_samplesize) <- c( "Control full", "Control sub", "UC full", "UC sub","CD full", "CD sub","IBD-U full")
fillscale_diagnosis_sample <- scale_fill_manual(name="Diagnosis",
                                                values = myColors_diagnosis_sample_samplesize, drop = T, limits=force)
colscale_diagnosis_sample <- scale_color_manual(name="Diagnosis",
                                                values = myColors_diagnosis_sample_samplesize, drop = T, limits=force)

