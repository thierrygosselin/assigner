## ------------------------------------------------------------------------
data = "batch_1.vcf" # or
data = "batch_1.haplotypes.tsv"
data = "tuna.tped" 
data = "tuna.tsv"

## ------------------------------------------------------------------------
mixture = "mixture.tuna.txt"

## ------------------------------------------------------------------------
whitelist.markers = "mixture.tuna.txt"

## ------------------------------------------------------------------------
blacklist.id = "blacklisted.individuals.txt"

## ------------------------------------------------------------------------
pop.levels = c("HAY", "ANG", "FEG")

## ------------------------------------------------------------------------
# to combined HAY and ANG into a new grouping: HAY-ANG
pop.levels = c("HAY", "ANG", "FEG")
pop.labels = c("HAY-ANG", "HAY-ANG", "FEG")

## ---- fig.show='hold'----------------------------------------------------
plot(1:10)
plot(10:1)

## ---- echo=FALSE, results='asis'-----------------------------------------
knitr::kable(head(mtcars, 10))

