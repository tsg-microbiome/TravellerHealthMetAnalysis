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

