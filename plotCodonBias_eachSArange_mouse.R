library(data.table)

readEachFile <- function(inputFile){
    return(read.csv(inputFile, header=F, col.names=c("POS","CODON","AA","SA","CN")))
}

getCodon2AAprobM <- function(allD, catNr){
    allD = allD[allD$SAcat==catNr,]
    return(prop.table(table(allD$AA, allD$CODON),1))
}

getCodon2AA_m <- function(wholeList, AA){
    noneZeroPos = which(wholeList[[1]][AA,] != 0)
    numCat = length(wholeList)

    tmpM = c()
    for (i in 1:numCat){
        tmpM = rbind(tmpM, wholeList[[i]][AA,][noneZeroPos])
    }
    rownames(tmpM) = 1:numCat

    return(tmpM)
}

#read in data
desDir="/home/kwang2/data/mouseData/MLRinput/inputs/CN_aa"
allDTs = lapply(Sys.glob(paste0(desDir, '/*.csv')), readEachFile)
allData = rbindlist(allDTs)

#determine SA category based on SA quantile (add a SAcat column)
#1                           [ 0.0,  0.3)
#2                           [ 0.3,  2.8)
#3                           [ 2.8,  8.6)
#4                           [ 8.6, 16.2)
#5                           [16.2, 25.5)
#6                           [25.5, 35.6)
#7                           [35.6, 45.8)
#8                           [45.8, 57.1)
#9                           [57.1, 72.0)
#10                           [72.0,172.2]
allData$SAcat = as.integer(cut2(allData$SA, g=10))

#in each SA category, tabulate codon given AA probs for each AA
codon2AAprobM_list = lapply(seq(max(allData$SAcat)), function(x) getCodon2AAprobM(allData, x))
#then for each AA, calculate a table like this:
#'A':
#  TTT TTA TTG
#  0.1 0.6 0.3 sum(...) = 1
#  ... ... ...
AAs = dimnames(codon2AAprobM_list[[1]])[[1]]

AAcodonProbM_list = lapply(AAs, function(x) getCodon2AA_m(codon2AAprobM_list,x))
names(AAcodonProbM_list) = AAs
#save AAcodonProbM_list so that I can try plot parameters on desktop
save(AAcodonProbM_list, file="AAcodonProbM_list.saved")

#plot part
for (i in 1:length(AAcodonProbM_list)){
    tmpAA = names(AAcodonProbM_list)[i]
    pdf(paste0(tmpAA, "_codons_SAcategories.pdf"))
    tmpM = AAcodonProbM_list[[i]]
    codonNum = ncol(tmpM)
    matplot(tmpM, type = c("b"),pch=1,col = 1:codonNum,
            ylab=paste0("Prob of AA type: ", tmpAA),
            xlab="SA category" )
    legend("topright", legend = 1:codonNum, col=1:codonNum, pch=1)
    dev.off()
}
