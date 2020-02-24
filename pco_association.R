#Species Comparison With Study as a confounder
library(lmtest)
library(QuantPsyc)

pco_association_with_study_confounder = function(data,pco,study_data)
{
	SpeciesAssociation <- as.data.frame(matrix(NA,ncol(data),2))
	rownames(SpeciesAssociation) <- colnames(data)
	for(i in 1:ncol(data))
	{
		print(i)
		fit1 <- lm(data[rownames(study_data),colnames(data)[i]]~as.factor(study_data$Study))
		fit2 <- lm(data[rownames(study_data),colnames(data)[i]]~as.factor(study_data$Study) + pco[rownames(study_data),2])
		t <- lrtest(fit1,fit2)
		SpeciesAssociation[i,2] <- t$Pr[2]
		fit3 <- lm(data[rownames(study_data),colnames(data)[i]]~pco[rownames(study_data),2])
		SpeciesAssociation[i,1] <- as.numeric(fit3$coefficients[2])
		i <- i + 1
	}
		return(SpeciesAssociation)
}
