# 2. 
# load the libraries
library(genomation)
library(GenomicRanges)
library(dplyr)
library(ggplot2)
library(ComplexHeatmap)

# load the data
chip.files = list.files('./Data/ChIP', full.names=TRUE)

lchip = lapply(chip.files, readGeneric, zero.based=TRUE)
lnams = basename(chip.files)
lnams = sub('_GRCh38.bed.gz','',lnams)
lnams = sub('K562_','',lnams)
names(lchip) = lnams
lchip = GRangesList(lchip)

# some samples have duplicated peaks and we need to collapse them
lchip = reduce(lchip)

# some samples have "werid" peak size distributions
uchip = unlist(lchip)
data.frame(sample = names(uchip), width = width(uchip)) %>%
    mutate(MAX = case_when(
        sample == 'MAX' ~ 'MAX',
        sample == 'TAL1' ~ 'TAL1',
        sample == 'FOXA1' ~ 'FOXA1',
        sample == 'ZNF639' ~ 'ZNF639',
        TRUE ~ 'Other'
        )) %>%
    ggplot(aes(width, color=MAX)) +
    geom_density(size=1) +
    scale_color_manual(values=c('darkorange','darkorange1','black','darkorange2','darkorange3','darkorange4')) +
    scale_x_log10() +
    theme_bw()

lchip = lchip[!names(lchip) %in% c('MAX','TAL1','FOXA1','ZNF639')]
uchip = unlist(lchip)

# intersection
ovlap      = as.data.frame(findOverlaps(uchip, lchip))
ovlap$set1 = names(uchip)[ovlap$queryHits]
ovlap$set2 = names(lchip)[ovlap$subjectHits]
inter      = ovlap %>%
    group_by(set1, set2) %>%
    summarize(intersection = length(unique(queryHits))) 
inter$key = with(inter, paste(set1, set2))

# union
union            = expand.grid(seq(lchip), seq(lchip))
union$set1       = names(lchip)[union$Var1]
union$set2       = names(lchip)[union$Var2]
union$length1    = elementNROWS(lchip)[union$Var1]
union$length2    = elementNROWS(lchip)[union$Var2]
union$sum_length = with(union, (length1 + length2))
union$key = with(union, paste(set1, set2))
union = union[,-c(1:4)]

# jaccard
dist = merge(inter, union, by='key')

# calculates the percentage of overlap for each TF - TF combination
dist$perc  = with(dist, round(intersection/length1,2))

# calculates the union
dist$union = with(dist, sum_length - intersection)

# calculates the approximation to the Jaccard index
dist$jacc  = with(dist, intersection/union)


mat.jacc = data.table::dcast(dist, set1~set2, value.var='jacc',fill=0) %>%
    mutate(set1 = NULL) 
rownames(mat.jacc) = colnames(mat.jacc)
mat.jacc = (mat.jacc + t(mat.jacc))/2


### Approximate Jaccard index
pdf('JaccardPeaks.pdf', width = 60, height=60)
    Heatmap(mat.jacc, col=circlize::colorRamp2(breaks = c(0,1),colors=c('white','red')))
dev.off()

### Average percentage of overlap

# heatmap
mat.perc = data.table::dcast(dist, set1~set2, value.var='perc',fill=0) %>%
    mutate(set1 = NULL) 
rownames(mat.perc) = colnames(mat.perc)
mat.perc = (mat.perc + t(mat.perc))/2


library(ComplexHeatmap)
pdf('PercentageHeatmap.pdf', width = 60, height=60)
    Heatmap(mat.perc, col=circlize::colorRamp2(breaks = c(0,1),colors=c('white','red')))
dev.off()

PairwiseJaccard = function(lchip, max_n = min(150, length(lchip))){
    
    uchip = unlist(lchip)
    d     = disjoin(uchip)
    dd    = data.frame(lapply(lchip[1:max_n], function(x)countOverlaps(d, x)>0))
    w     = width(d)
    ls = lapply(head(seq(dd),-1), function(x){
        lapply(seq(x+1, ncol(dd)), function(y){
            message(paste(x,y))
            sum(w[dd[,x] & dd[,y]])
        })
    })
    ls = lapply(ls, unlist)
    m = matrix(0,ncol=max_n, nrow=max_n)
    m[lower.tri(m)] = unlist(ls)
    m = m + t(m)
    wl = sum(width(lchip[1:max_n]))
    o = outer(wl, wl, '+') - m

    jacc = m/o
    diag(jacc) = 1
    round(jacc,2)
}
gr  = GRanges('chr1', IRanges(c(1,10), width=4))

toy = GRangesList(gr, shift(gr,2), shift(gr, 4), shift(gr, 5))

jacc = PairwiseJaccard(toy)


pdf('Jaccard_Heatmap.pdf', width = 30, height=30)
    Heatmap(jacc, col=circlize::colorRamp2(breaks = c(0,1),colors=c('white','red')))
dev.off()

