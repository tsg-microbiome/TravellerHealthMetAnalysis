#PERMANOVA Global
#Objects: NewCombinedSpecies, dfStudyPopulationTypeGlobal, SchirmerM_2016Individuals, Industrialized, NonIndustrialized, NonIndustrializedLike, Intermediate, IndustrializedLike, all_other_irish

AdonisGlobalWithStudyConfounding <- adonis(as.dist(1-cor(t(NewCombinedSpecies[setdiff(rownames(dfStudyPopulationTypeGlobal),SchirmerM_2016Individuals),]),method="spearman")/2)~dfStudyPopulationTypeGlobal[setdiff(rownames(dfStudyPopulationTypeGlobal),SchirmerM_2016Individuals),"Study"]+dfStudyPopulationTypeGlobal[setdiff(rownames(dfStudyPopulationTypeGlobal),SchirmerM_2016Individuals),"Population"])
lm_pcoAxis1_WithPopulation <- anova(lm(pcoNewCombinedSpecies$li[setdiff(rownames(dfStudyPopulationTypeGlobal),SchirmerM_2016Individuals),1]~dfStudyPopulationTypeGlobal[setdiff(rownames(dfStudyPopulationTypeGlobal),SchirmerM_2016Individuals),"Study"]+dfStudyPopulationTypeGlobal[setdiff(rownames(dfStudyPopulationTypeGlobal),SchirmerM_2016Individuals),"Population"]))
lm_pcoAxis2_WithPopulation <- anova(lm(pcoNewCombinedSpecies$li[setdiff(rownames(dfStudyPopulationTypeGlobal),SchirmerM_2016Individuals),2]~dfStudyPopulationTypeGlobal[setdiff(rownames(dfStudyPopulationTypeGlobal),SchirmerM_2016Individuals),"Study"]+dfStudyPopulationTypeGlobal[setdiff(rownames(dfStudyPopulationTypeGlobal),SchirmerM_2016Individuals),"Population"]))
lm_pcoAxis2_WithPopulation <- anova(lm(pcoNewCombinedSpecies$li[setdiff(rownames(dfStudyPopulationTypeGlobal),SchirmerM_2016Individuals),3]~dfStudyPopulationTypeGlobal[setdiff(rownames(dfStudyPopulationTypeGlobal),SchirmerM_2016Individuals),"Study"]+dfStudyPopulationTypeGlobal[setdiff(rownames(dfStudyPopulationTypeGlobal),SchirmerM_2016Individuals),"Population"]))

#PERMANOVA Irish
#Objects: IrishSampleTags, IrishSampleTagsWithStudy
AdonisIrishWithStudyConfounding <- adonis(as.dist(1-cor(t(NewCombinedSpecies[rownames(IrishSampleTagsWithStudy),]),method="spearman")/2)~IrishSampleTagsWithStudy[,"Study"]+IrishSampleTagsWithStudy[,"IrishSampleTags"])
AdonisIrish <- adonis(as.dist(1-cor(t(NewCombinedSpecies[rownames(IrishSampleTagsWithStudy),]),method="spearman")/2)~IrishSampleTagsWithStudy[,"IrishSampleTags"])

#Mann-Whitney Comparison Across Groups
batch_wilcox <- function(x,y)
{
        p_array <- NULL;
        type_array <- NULL;
        mean1_array <- NULL;
        mean2_array <- NULL;
        #x <- x[abs(rowSums(x,na.rm=TRUE)) > 0,];
        #y <- y[abs(rowSums(y,na.rm=TRUE)) > 0,];
        z <- intersect(rownames(x),rownames(y));
        for(i in 1:length(z))
        {
                p_array[i] <- wilcox.test(as.numeric(x[z[i],]),as.numeric(y[z[i],]))$p.value;
                type_array[i] <- ifelse(mean(as.numeric(x[z[i],]),na.rm=TRUE) > mean(as.numeric(y[z[i],]),na.rm=TRUE), 1, ifelse(mean(as.numeric(x[z[i],]),na.rm=TRUE) < mean(as.numeric(y[z[i],]),na.rm=TRUE),-1,0));
                mean1_array[i] <- mean(as.numeric(x[z[i],]),na.rm=TRUE);
                mean2_array[i] <- mean(as.numeric(y[z[i],]),na.rm=TRUE);
                i <- i + 1;
        }
        out <- as.data.frame(cbind(p_array,type_array,p.adjust(p_array),mean1_array,mean2_array));
        rownames(out) <- z;
        out <- apply(out,1,function(x)(ifelse(is.nan(x),1,x)));
        return(t(out));
}

SelectSpecies <- colnames(NewCombinedSpecies)[(which(100*apply(NewCombinedSpecies,2,function(x)(length(x[x>=0.01])))/ncol(NewCombinedSpecies)>5))]
BatchWilcoxSpeciesNonIndustrialized_Industrialized <- batch_wilcox(t(NewCombinedSpecies[setdiff(Industrialized,SchirmerM_2016Individuals),SelectSpecies]),t(NewCombinedSpecies[NonIndustrialized,SelectSpecies]))
SpeciesEnrichedNonIndustrialized <- rownames(BatchWilcoxSpeciesNonIndustrialized_Industrialized[(BatchWilcoxSpeciesNonIndustrialized_Industrialized[,2]==-1)&(BatchWilcoxSpeciesNonIndustrialized_Industrialized[,3]<0.15),])
SpeciesDepletedNonIndustrialized <- rownames(BatchWilcoxSpeciesNonIndustrialized_Industrialized[(BatchWilcoxSpeciesNonIndustrialized_Industrialized[,2]==1)&(BatchWilcoxSpeciesNonIndustrialized_Industrialized[,3]<0.15),])

SelectSpecies_TM <- colnames(tm_species)[which(100*apply(tm_species,2,function(x)(length(x[x>=0.01])))/ncol(tm_species)>5)]
BatchWilcoxSpeciesNonIndustrializedLike_IndustrializedLike <- batch_wilcox(t(NewCombinedSpecies[IndustrializedLike,SelectSpecies_TM]),t(NewCombinedSpecies[NonIndustrializedLike,SelectSpecies_TM]))
SpeciesEnrichedNonIndustrializedLike <- rownames(BatchWilcoxSpeciesNonIndustrializedLike_IndustrializedLike[(BatchWilcoxSpeciesNonIndustrializedLike_IndustrializedLike[,2]==-1)&(BatchWilcoxSpeciesNonIndustrializedLike_IndustrializedLike[,3]<0.15),])
SpeciesDepletedNonIndustrializedLike <- rownames(BatchWilcoxSpeciesNonIndustrializedLike_IndustrializedLike[(BatchWilcoxSpeciesNonIndustrializedLike_IndustrializedLike[,2]==1)&(BatchWilcoxSpeciesNonIndustrializedLike_IndustrializedLike[,3]<0.15),])


#Species Validation With Study as a confounder
library(lmtest)
library(QuantPsyc)

pco_association_with_study_confounder = function(data,pco,study_data)
{
	SpeciesAssociation <- as.data.frame(matrix(NA,ncol(data),2))
	rownames(SpeciesAssociation) <- colnames(data)
	for(i in 1:ncol(data))
	{
		fit1 <- lm(data[rownames(study_data),colnames(data)[i]]~study_data$Study)
		fit2 <- lm(data[rownames(study_data),colnames(data)[i]]~study_data$Study + pco[rownames(study_data),2])
		t <- lrtest(fit1,fit2)
		SpeciesAssociation[i,2] <- t$Pr[2]
		fit3 <- lm(data[rownames(study_data),colnames(data)[i]]~pco[rownames(study_data),2])
		SpeciesAssociation[i,1] <- as.numeric(fit3$coefficients[2])
		i <- i + 1
	}
		return(SpeciesAssociation)
}

pco2AssociationSpeciesDepleted <- pco_association_with_study_confounder(NewCombinedSpecies[setdiff(c(NonIndustrialized,Industrialized),SchirmerM_2016Individuals),SpeciesDepletedNonIndustrialized],pcoNewCombinedSpecies$li[setdiff(c(NonIndustrialized,Industrialized),SchirmerM_2016Individuals),],dfStudyPopulationTypeGlobal[setdiff(c(NonIndustrialized,Industrialized),SchirmerM_2016Individuals),])
pco2AssociationSpeciesEnriched <- pco_association_with_study_confounder(NewCombinedSpecies[setdiff(c(NonIndustrialized,Industrialized),SchirmerM_2016Individuals),SpeciesEnrichedNonIndustrialized],pcoNewCombinedSpecies$li[setdiff(c(NonIndustrialized,Industrialized),SchirmerM_2016Individuals),],dfStudyPopulationTypeGlobal[setdiff(c(NonIndustrialized,Industrialized),SchirmerM_2016Individuals),])

#Distance Profile of Microbiome Sub-groups and Residence Transition Groups 
#Object: MicrobiomeResidenceTransitions Site_Site Site_House House_House AllIrish_NIvsI_Distance_Profile
library(dunn.test)

dunn.test(MicrobiomeResidenceTransitions[,"NonIndustrialized"],as.factor(MicrobiomeResidenceTransitions[,3],method="bh")
dunn.test(MicrobiomeResidenceTransitions[,"Industrialized"],as.factor(MicrobiomeResidenceTransitions[,3],method="bh")
dunn.test(AllIrish_NIvsI_Distance_Profile[,1],as.factor(AllIrish_NIvsI_Distance_Profile[,3]),method="bh")
dunn.test(AllIrish_NIvsI_Distance_Profile[,2],as.factor(AllIrish_NIvsI_Distance_Profile[,3]),method="bh")

#Objects: CombinedHeFD
#HeFD Comparisons
CombinedHeFDTags <- factor(c(rep("NonIndustrializedLike",length(NonIndustrializedLike)),rep("Intermediate",length(Intermediate)),rep("IndustrializedLike",length(IndustrializedLike)),rep("EMLongstay",length(EMLongstay)),rep("EMCommunity",length(EMCommunity)),rep("YoungIrish",length(GMHeFD))),levels=c("NonIndustrializedLike","Intermediate","IndustrializedLike","EMLongstay","EMCommunity","YoungIrish"))
dunn.test(CombinedHeFD,as.factor(CombinedHeFDTags))

#TMA Producer Analysis
#Objects: TMA_Producers, TMA_Consumers, Final_CorrP_TMAProducers_WithFood, CombinedFoodCodeGrouped
TMA_Producers_WithinTravellers <- batch_wilcox(t(NewCombinedSpecies[NonIndustrializedLike,intersect(TMA_Producers,colnames(NewCombinedSpecies))]),t(NewCombinedSpecies[IndustrializedLike,intersect(TMA_Producers,colnames(NewCombinedSpecies))]))
Detected_TMA_Producers <- rownames(TMA_Producers_WithinTravellers[TMA_Producers_WithinTravellers[,2]!=0,])
heatmap.2(t(Final_CorrP_TMAProducers_WithFood$r[Detected_TMA_Producers,]),density="none",trace="none",col=brewer.pal(8,"PiYG"),cellnote=apply(Final_CorrP_TMAProducers_WithFood$p,1,function(x)(ifelse(x<0.05,"*",""))),notecol="black",notecex=1.5,margins=c(15,25),lhei=c(1,5),lwid=c(1,5),cexRow=1.2)

CorrP_NI_I_MajorFoodGroup <- corr.p(cor(AllMetaTable[,c("NonIndustrialized","Industrialized")],CombinedMacroDietNormalized[rownames(AllMetaTable),SelectFoodCategories],method="spearman"),n=117,adjust="fdr")
heatmap.2(CorrP_NI_I_MajorFoodGroup$r,density="none",trace="none",col=brewer.pal(8,"RdYlBu"),margins=c(35,25),Rowv=FALSE,lhei=c(1,5),lwid=c(1,5),cexCol=4,cellnote=apply(CorrP_NI_I_MajorFoodGroup$p,2,function(x)(ifelse(x<0.1,"*",""))),notecex=3,notecol="black")

#Clinical Metadata Analysis
#Objects: CorrP_NonIndustrialized, CorrP_Industrialized
dfClinical_NonIndustrialized <- as.data.frame(cbind(CorrP_NonIndustrialized$r,CorrP_NonIndustrialized$p))
colnames(dfClinical_NonIndustrialized) <- c("R","P")
dfClinical_Industrialized <- as.data.frame(cbind(CorrP_Industrialized$r,CorrP_Industrialized$p))
colnames(dfClinical_Industrialized) <- c("R","P")
dfClinical_NonIndustrialized <- dfClinical_NonIndustrialized[setdiff(rownames(dfClinical_NonIndustrialized),"Creatinine"),]
dfClinical_Industrialized <- dfClinical_Industrialized[setdiff(rownames(dfClinical_Industrialized),"Creatinine"),]
ggplot(dfClinical_NonIndustrialized,aes(x=R,y=-log(P,10))) + geom_point(color=ifelse(dfClinical_NonIndustrialized$P < 0.1,"Red",ifelse(dfClinical_NonIndustrialized$P < 0.2,"Darkgoldenrod1","Black")),size=10) + geom_text_repel(label=rownames(dfClinical_NonIndustrialized),size=10) + theme_bw() + geom_vline(xintercept=0) + geom_hline(yintercept=0) + theme(axis.text=element_text(size=20))
ggplot(dfClinical_Industrialized,aes(x=R,y=-log(P,10))) + geom_point(color=ifelse(dfClinical_Industrialized$P < 0.1,"Red",ifelse(dfClinical_Industrialized$P < 0.2,"Darkgoldenrod1","Black")),size=10) + geom_text_repel(label=rownames(dfClinical_Industrialized),size=10) + theme_bw() + geom_vline(xintercept=0) + geom_hline(yintercept=0) + theme(axis.text=element_text(size=20))

#Objects: 

