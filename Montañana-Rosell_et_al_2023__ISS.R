
##############################################
#### In situ sequencing analysis pipeline ####
##############################################

### Prep R

library(readxl)
library(dplyr)
library(tidyr)
library(gdata)
library(rearrr)
library(ggplot2)
library(ggthemes)



#######################
####   PREP DATA   ####
#######################

### Load reads files

P1_1 <- read.csv("P1_LUMBAR_4-1_0_4.csv")
P1_2 <- read.csv("P1_LUMBAR_4-2_0_4.csv")
P1_3 <- read.csv("P1_LUMBAR_4-3_0_4.csv")
P1_4 <- read.csv("P1_LUMBAR_4-4_0_4.csv")
P1_5 <- read.csv("P1_LUMBAR_4-5_0_4.csv")
P1_6 <- read.csv("P1_LUMBAR_4-6_0_4.csv")
P28_1 <- read.csv("section1_Z3_0_4.csv")
P28_2 <- read.csv("section2_Z4_0_4.csv")
P28_3 <- read.csv("P28_Gama_7-1_0_4.csv")
P28_4 <- read.csv("P28_Gama_7-2_0_4.csv")
P28_5 <- read.csv("P28_Gama_7-3_0_4.csv")
P28_C6.1 <- read.csv("GASTRO_C6-1_0_4.csv")
P28_C6.2 <- read.csv("GASTRO_C6-2_0_4.csv")
P28_C6.3 <- read.csv("GASTRO_C6-3_0_4.csv")
P28_C6.4 <- read.csv("GASTRO_C6-4_0_4.csv")
P28_C6.5 <- read.csv("GASTRO_C6-5_0_4.csv")
P28_C6.6 <- read.csv("GASTRO_C6-6_0_4.csv")
P28_C7.1 <- read.csv("GASTRO_C7-1_0_4.csv")
P28_C7.2 <- read.csv("GASTRO_C7-2_0_4.csv")
P28_C7.3 <- read.csv("GASTRO_C7-3_0_4.csv")
P28_C7.4 <- read.csv("GASTRO_C7-4_0_4.csv")
P28_C7.5 <- read.csv("GASTRO_C7-5_0_4.csv")
P28_C7.6 <- read.csv("GASTRO_C7-6_0_4.csv")
P28_E3.1 <- read.csv("SOLEUS_E3-1_0_4.csv")
P28_E3.2 <- read.csv("SOLEUS_E3-2_0_4.csv")
P28_E3.3 <- read.csv("SOLEUS_E3-3_0_4.csv")
P28_E3.4 <- read.csv("SOLEUS_E3-4_0_4.csv")
P28_E3.5 <- read.csv("SOLEUS_E3-5_0_4.csv")
P28_E3.6 <- read.csv("SOLEUS_E3-6_0_4.csv")
P28_E4.1 <- read.csv("SOLEUS_E4-1_0_4.csv")
P28_E4.2 <- read.csv("SOLEUS_E4-2_0_4.csv")
P28_E4.3 <- read.csv("SOLEUS_E4-3_0_4.csv")
P28_E4.4 <- read.csv("SOLEUS_E4-4_0_4.csv")
P28_E4.5 <- read.csv("SOLEUS_E4-5_0_4.csv")
P28_E4.6 <- read.csv("SOLEUS_E4-6_0_4.csv")


### Remove Chx10 (Vsx2) for 2021 samples due to unspecific signal

P28_C6.1 <- subset(P28_C6.1, gene!="Vsx2")
P28_C6.2 <- subset(P28_C6.2, gene!="Vsx2")
P28_C6.3 <- subset(P28_C6.3, gene!="Vsx2")
P28_C6.4 <- subset(P28_C6.4, gene!="Vsx2")
P28_C6.5 <- subset(P28_C6.5, gene!="Vsx2")
P28_C6.6 <- subset(P28_C6.6, gene!="Vsx2")
P28_C7.1 <- subset(P28_C7.1, gene!="Vsx2")
P28_C7.2 <- subset(P28_C7.2, gene!="Vsx2")
P28_C7.3 <- subset(P28_C7.3, gene!="Vsx2")
P28_C7.4 <- subset(P28_C7.4, gene!="Vsx2")
P28_C7.5 <- subset(P28_C7.5, gene!="Vsx2")
P28_C7.6 <- subset(P28_C7.6, gene!="Vsx2")
P28_E3.1 <- subset(P28_E3.1, gene!="Vsx2")
P28_E3.2 <- subset(P28_E3.2, gene!="Vsx2")
P28_E3.3 <- subset(P28_E3.3, gene!="Vsx2")
P28_E3.4 <- subset(P28_E3.4, gene!="Vsx2")
P28_E3.5 <- subset(P28_E3.5, gene!="Vsx2")
P28_E3.6 <- subset(P28_E3.6, gene!="Vsx2")
P28_E4.1 <- subset(P28_E4.1, gene!="Vsx2")
P28_E4.2 <- subset(P28_E4.2, gene!="Vsx2")
P28_E4.3 <- subset(P28_E4.3, gene!="Vsx2")
P28_E4.4 <- subset(P28_E4.4, gene!="Vsx2")
P28_E4.5 <- subset(P28_E4.5, gene!="Vsx2")
P28_E4.6 <- subset(P28_E4.6, gene!="Vsx2")



### Merge all datasets (sample of origin will be named as "source") and plot by sample to preview

ISS <- combine(P1_1, P1_2, P1_3, P1_4, P1_5, P1_6,
               P28_1, P28_2, P28_3, P28_4, P28_5,
               P28_C6.1, P28_C6.2, P28_C6.3, P28_C6.4, P28_C6.5, P28_C6.6,
               P28_C7.1, P28_C7.2, P28_C7.3, P28_C7.4, P28_C7.5, P28_C7.6,
               P28_E3.1, P28_E3.2, P28_E3.3, P28_E3.4, P28_E3.5, P28_E3.6,
               P28_E4.1, P28_E4.2, P28_E4.3, P28_E4.4, P28_E4.5, P28_E4.6)

# ISS <- combine(P1_1, P1_2, P28_1, P28_2)

ggplot(ISS, aes(X, Y)) +
  geom_point(size = 1) +
  facet_wrap(~ source) +
  theme_few() +
  ggtitle ("Original coordinates")




###########################
####   NORMALIZATION   ####
###########################


### ROTATION

## Load .xlsx file for CC coordinates and angle
## Add rotation info to dataframe based on source
## Transform X, Y coordinates to new 0,0 by subtraction of CC coordinates (and delete CC coordinates)
## Rotate coordinates based on angle and new 0,0 origin using group_map to split the data by source and separately apply the rotate_2d function, then bind and delete extra columns
## (for the angle, use either mean/unique/min to pick a single value - doesn't matter since they are all the same)
## Visualize cells after rotation

Rotation <- read_xlsx("Rotation&CC.xlsx")

ISS <- merge(ISS, Rotation, by="source", all.x=TRUE)

ISS$X <- (ISS$X - ISS$CC_X)
ISS$Y <- (ISS$Y - ISS$CC_Y)
ISS <- subset(ISS, select = -c(CC_X, CC_Y))

ISS <- ISS %>% group_by(source) %>% 
  group_map(
    ~ rotate_2d(data = ., degrees = mean(.[["Angle"]], na.rm = TRUE), x_col = "X", y_col = "Y", origin = c(0,0), suffix = '', overwrite = TRUE),
    .keep = TRUE
  ) %>% bind_rows()

ISS <- subset(ISS, select = -c(Angle, .origin, .degrees))
rm(Rotation)

ggplot(ISS, aes(X, Y)) +
  geom_point(size = 1) +
  facet_wrap(~ source) +
  theme_few() +
  ggtitle("Coordinates after rotation")



### SPLIT LEFT/RIGHT

## Add column stating left/right depending on X coordinate, change X vlues to positives for left hemicords, rename "source" to include hemicord info
## Visualize cells after hemicord split

ISS %>% mutate(Hemicord = case_when(X >= 0 ~ "Right", X < 0 ~ "Left")) -> ISS

ISS %>% mutate(X = case_when(Hemicord == "Left" ~ X*(-1), Hemicord == "Right" ~ X)) -> ISS

ISS %>% unite(source, c(source, Hemicord), sep = "_") -> ISS

ggplot(ISS, aes(X, Y)) +
  geom_point(size = 1) +
  facet_wrap(~ source, nrow = 6) +
  theme_few() +
  ggtitle("Coordinates after hemicord split")



### NORMALIZE FOR WIDTH/HEIGHT

## Load .xlsx file for width and height normalization
## Add columns for width & height factors and normalize X,Y coordinates, then visualize

Size <- read_xlsx("Size_normalization.xlsx")

ISS <- merge(ISS, Size %>% select(source, Width, Height), by="source", all.x=TRUE)

ISS$X <- (ISS$X * ISS$Width)
ISS$Y <- (ISS$Y * ISS$Height)

ISS <- subset(ISS, select = -c(Width, Height))

rm(Size)

ggplot(ISS, aes(X, Y)) +
  geom_point(size = 1) +
  facet_wrap(~ source, nrow = 6) +
  theme_few() +
  ggtitle("Coordinates after width/height normalization")



### CHANGE FROM pixel TO um (0.32um/pixel)

ISS$X <- (ISS$X * 0.32)
ISS$Y <- (ISS$Y * 0.32)



### SPLIT SOURCE INFO TO TIMEPOINT/SAMPLE/HEMICORD (and visualize by gene)

ISS <- separate(ISS, source, c("Timepoint", "source", "L/R"), "_")
ISS <- ISS %>% unite(Hemisection, c("source", "L/R"), sep = "_")

ggplot(ISS, aes(X, Y)) +
  xlim(0, 1500) + xlab("X (um)") +
  ylab("Y (um)") + 
  geom_point(aes(color = gene), size = 1) +
  facet_wrap(Timepoint ~ .) +
  guides(color = guide_legend(ncol = 1)) +
  theme_few() +
  ggtitle("Coordinates after normalization")




######################
####   PLOTTING   ####
######################


### RENAME GENES TO MORE COMMON NAMES

ISS[ISS == "Bhlhe22"] <- "Bhlhb5"
# ISS[ISS == "Calb1"] <- "Calbindin"
# ISS[ISS == "Calb2"] <- "Calretinin"
ISS[ISS == "Chat"] <- "ChAT"
ISS[ISS == "Gad1"] <- "Gad67"
ISS[ISS == "Gad2"] <- "Gad65"
ISS[ISS == "Slc17a6"] <- "Vglut2"
ISS[ISS == "Slc6a5"] <- "GlyT2"
ISS[ISS == "Vsx2"] <- "Chx10"
ISS <- with(ISS,  ISS[order(gene) , ])



### PLOT
## All samples, all markers

windowsFonts(Arial = windowsFont("Arial"))

ggplot(ISS, aes(X, Y)) +
  scale_x_continuous(limits = c(0, 1100), breaks = seq(0, 1000, 500)) + xlab("X (um)") +
  scale_y_continuous(limits = c(-700, 700), breaks = seq(-500, 500, 500)) + ylab("Y (um)") + 
  geom_point(aes(color = gene), size = 1) +
  facet_wrap(Timepoint ~ .) +
  guides(color = guide_legend(ncol = 1)) +
  theme_few(base_family = "Arial") +
  theme(axis.title.y = element_text(color = "black", size = 14), axis.text.y = element_text(color = "black", size = 14),
        axis.title.x = element_text(color = "black", size = 14), axis.text.x = element_text(color = "black", size = 14),
        axis.ticks.y = element_line(color = "black"), axis.ticks.x = element_line(color = "black"),
        panel.border = element_rect(color = "black", size = 1), strip.text.x = element_text(color = "black", size = 18))


## Select specific subset of data to display
#(e.g. L&R hemisections for first two samples of P1 & P28 -at this point they will have the same Hemisection name-)

# Samples 1 & 2 for P1 & P28
ISS_fig1 <- subset(ISS, Hemisection == "1_Left" | Hemisection == "1_Right" | Hemisection == "2_Left" | Hemisection == "2_Right")
# Samples 1-6 for P1, 1-5 for P28
ISS_fig2 <- subset(ISS, Hemisection == "1_Left" | Hemisection == "1_Right" | Hemisection == "2_Left" | Hemisection == "2_Right" |
                     Hemisection == "3_Left" | Hemisection == "3_Right" | Hemisection == "4_Left" | Hemisection == "4_Right" |
                     Hemisection == "5_Left" | Hemisection == "5_Right" | Hemisection == "6_Left" | Hemisection == "6_Right")


## Select markers of interest

ggplot(filter(ISS_fig1, gene == "En1" | gene == "Chx10" | gene == "Shox2" | gene == "Pitx2" | gene == "ChAT")) +
#ggplot(filter(ISS_fig1, gene == "Calb1" | gene == "Calb2" | gene == "Foxp2" | gene == "Pou6f2" | gene == "Sp8")) +
  scale_x_continuous(limits = c(0, 1100), breaks = seq(0, 1000, 500)) + xlab("X (um)") +
  scale_y_continuous(limits = c(-700, 700), breaks = seq(-500, 500, 500)) + ylab("Y (um)") + 
  geom_point(aes(X, Y, color = gene), size = 2) +
  facet_wrap(Timepoint ~ .) +
  theme_few(base_family = "Arial") +
  theme(axis.title.y = element_text(color = "black", size = 14), axis.text.y = element_text(color = "black", size = 14),
        axis.title.x = element_text(color = "black", size = 14), axis.text.x = element_text(color = "black", size = 14),
        axis.ticks.y = element_line(color = "black"), axis.ticks.x = element_line(color = "black"),
        panel.border = element_rect(color = "black", size = 1), strip.text.x = element_text(color = "black", size = 18))




#########################
#### CO-LOCALIZATION ####
#########################

### OPTIONAL, TO UNDERSTAND HOW MANY CO-LOCALIZATIONS THERE ARE AND HOW MANY CELLS HAVE ONLY ONE TRANSCRIPT
## Subset only cells that appear more than once (duplicated), or only once (unique)
#library(misty)
#ISS_duplicated <- df.duplicated(ISS, "X", "Y")
#ISS_unique <- df.unique(ISS, "Hemisection", "Timepoint")



### PREP: SPLIT DATAFRAME BY MARKERSplit dataframe by marker

ISS_marker <- split(ISS, ISS$gene)                ## All samples
#ISS_marker <- split(ISS_fig1, ISS_fig1$gene)     ## 2 samples
#ISS_marker <- split(ISS_fig2, ISS_fig2$gene)     ## 6 P1 samples + 5 P28 samples



### CREATE DATASETS FOR EACH CO-LOCALIZATION
## V1 co-localizations

En1_Foxp2 <- merge(ISS_marker$En1, ISS_marker$Foxp2, by = c("Timepoint", "Hemisection", "Area", "X", "Y"), all = FALSE)
En1_Pou6f2 <- merge(ISS_marker$En1, ISS_marker$Pou6f2, by = c("Timepoint", "Hemisection", "Area", "X", "Y"), all = FALSE)
En1_Sp8 <- merge(ISS_marker$En1, ISS_marker$Sp8, by = c("Timepoint", "Hemisection", "Area", "X", "Y"), all = FALSE)
En1_Calb1 <- merge(ISS_marker$En1, ISS_marker$Calb1, by = c("Timepoint", "Hemisection", "Area", "X", "Y"), all = FALSE)
En1_Calb2 <- merge(ISS_marker$En1, ISS_marker$Calb2, by = c("Timepoint", "Hemisection", "Area", "X", "Y"), all = FALSE)
#En1_Pvalb <- merge(ISS_marker$En1, ISS_marker$Pvalb, by = c("Timepoint", "Hemisection", "Area", "X", "Y"), all = FALSE)

V1 <- rbind(En1_Foxp2, En1_Pou6f2, En1_Sp8, En1_Calb1, En1_Calb2)
V1 <- V1 %>% unite(gene, c("gene.x", "gene.y"), sep = "/")

V1 <- rbind(ISS_marker$En1, V1)

ggplot(V1, aes(X, Y)) +
  scale_x_continuous(limits = c(0, 1100), breaks = seq(0, 1000, 500)) + xlab("X (um)") +
  scale_y_continuous(limits = c(-700, 700), breaks = seq(-500, 500, 500)) + ylab("Y (um)") + 
  geom_point(aes(color = gene), size = 2) +
  facet_wrap(Timepoint ~ .) +
  guides(color = guide_legend(ncol = 1)) +
  theme_few(base_family = "Arial") +
  theme(axis.title.y = element_text(color = "black", size = 14), axis.text.y = element_text(color = "black", size = 14),
        axis.title.x = element_text(color = "black", size = 14), axis.text.x = element_text(color = "black", size = 14),
        axis.ticks.y = element_line(color = "black"), axis.ticks.x = element_line(color = "black"),
        panel.border = element_rect(color = "black", size = 1), strip.text.x = element_text(color = "black", size = 18))



## V2a/Shox2 co-localizations

Chx10_Shox2 <- merge(ISS_marker$Chx10, ISS_marker$Shox2, by = c("Timepoint", "Hemisection", "Area", "X", "Y"), all = FALSE)
V2a <- Chx10_Shox2 %>% unite(gene, c("gene.y", "gene.x"), sep = "/")
V2a <- rbind(ISS_marker$Chx10, ISS_marker$Shox2, V2a)

ggplot(V2a, aes(X, Y)) +
  scale_x_continuous(limits = c(0, 1100), breaks = seq(0, 1000, 500)) + xlab("X (um)") +
  scale_y_continuous(limits = c(-700, 700), breaks = seq(-500, 500, 500)) + ylab("Y (um)") + 
  geom_point(aes(color = gene), size = 2) +
  facet_wrap(Timepoint ~ .) +
  guides(color = guide_legend(ncol = 1)) +
  theme_few(base_family = "Arial") +
  theme(axis.title.y = element_text(color = "black", size = 14), axis.text.y = element_text(color = "black", size = 14),
        axis.title.x = element_text(color = "black", size = 14), axis.text.x = element_text(color = "black", size = 14),
        axis.ticks.y = element_line(color = "black"), axis.ticks.x = element_line(color = "black"),
        panel.border = element_rect(color = "black", size = 1), strip.text.x = element_text(color = "black", size = 18))



## Inhibitory neurotransmitter co-localizations

GlyT2_Gad67 <- merge(ISS_marker$GlyT2, ISS_marker$Gad67, by = c("Timepoint", "Hemisection", "Area", "X", "Y"), all = FALSE)
GlyT2_Gad65 <- merge(ISS_marker$GlyT2, ISS_marker$Gad65, by = c("Timepoint", "Hemisection", "Area", "X", "Y"), all = FALSE)
Gad67_Gad65 <- merge(ISS_marker$Gad67, ISS_marker$Gad65, by = c("Timepoint", "Hemisection", "Area", "X", "Y"), all = FALSE)

Inhibitory <- rbind(GlyT2_Gad67, GlyT2_Gad65, Gad67_Gad65)
Inhibitory <- Inhibitory %>% unite(gene, c("gene.x", "gene.y"), sep = "/")

GlyT2_Gad67_Gad65 <- merge(GlyT2_Gad67, ISS_marker$Gad65, by = c("Timepoint", "Hemisection", "Area", "X", "Y"), all = FALSE)
GlyT2_Gad67_Gad65 <- GlyT2_Gad67_Gad65 %>% unite(gene, c("gene.x", "gene.y", "gene"), sep = "/")

Inhibitory <- rbind(ISS_marker$GlyT2, ISS_marker$Gad67, ISS_marker$Gad65, Inhibitory, GlyT2_Gad67_Gad65)

Inhibitory$gene <- factor(Inhibitory$gene, levels = c("GlyT2", "Gad67", "Gad65", "Gad67/Gad65", "GlyT2/Gad67", "GlyT2/Gad65", "GlyT2/Gad67/Gad65"))

ggplot(Inhibitory, aes(X, Y)) +
  scale_x_continuous(limits = c(0, 1100), breaks = seq(0, 1000, 500)) + xlab("X (um)") +
  scale_y_continuous(limits = c(-700, 700), breaks = seq(-500, 500, 500)) + ylab("Y (um)") + 
  geom_point(aes(color = gene), size = 2) +
  facet_wrap(Timepoint ~ .) +
  guides(color = guide_legend(ncol = 1)) +
  theme_few(base_family = "Arial") +
  theme(axis.title.y = element_text(color = "black", size = 14), axis.text.y = element_text(color = "black", size = 14),
        axis.title.x = element_text(color = "black", size = 14), axis.text.x = element_text(color = "black", size = 14),
        axis.ticks.y = element_line(color = "black"), axis.ticks.x = element_line(color = "black"),
        panel.border = element_rect(color = "black", size = 1), strip.text.x = element_text(color = "black", size = 18))



## Excitatory neurotransmitter co-localizations

Vglut2_ChAT <- merge(ISS_marker$Vglut2, ISS_marker$ChAT, by = c("Timepoint", "Hemisection", "Area", "X", "Y"), all = FALSE)
Excitatory <- Vglut2_ChAT %>% unite(gene, c("gene.x", "gene.y"), sep = "/")

Excitatory <- rbind(ISS_marker$Vglut2, ISS_marker$ChAT, Excitatory)

Excitatory$gene <- factor(Excitatory$gene, levels = c("Vglut2", "ChAT", "Vglut2/ChAT"))

ggplot(Excitatory, aes(X, Y)) +
  scale_x_continuous(limits = c(0, 1100), breaks = seq(0, 1000, 500)) + xlab("X (um)") +
  scale_y_continuous(limits = c(-700, 700), breaks = seq(-500, 500, 500)) + ylab("Y (um)") + 
  geom_point(aes(color = gene), size = 2) +
  facet_wrap(Timepoint ~ .) +
  guides(color = guide_legend(ncol = 1)) +
  theme_few(base_family = "Arial") +
  theme(axis.title.y = element_text(color = "black", size = 14), axis.text.y = element_text(color = "black", size = 14),
        axis.title.x = element_text(color = "black", size = 14), axis.text.x = element_text(color = "black", size = 14),
        axis.ticks.y = element_line(color = "black"), axis.ticks.x = element_line(color = "black"),
        panel.border = element_rect(color = "black", size = 1), strip.text.x = element_text(color = "black", size = 18))


