

# Import CBNA relevés CCVCMB   -------------------------------------------

data <- readRDS("_data/cbna_phyto_chamonix.rds")

sumary(data)
#les reléves datent que de 2017-2018 ?

length(unique(data$numchrono)) #seulement 80 relevés !!



# Identify distribution of key LANDE sp. --------------------------------------------------


#Rhododendron ferrugineum L.    117679
#Vaccinium myrtillus L.         128345
#Vaccinium uliginosum L.        128354
#Juniperus sibirica             104414  (à vérifier avec J. communis nana)

rhodo <- data[which(data$taxref_cd_ == 117679),]
hist(rhodo$rel_surfac)

myrt <- data[which(data$taxref_cd_ == 128345),]
hist(myrt$rel_surfac)

ulig <- data[which(data$taxref_cd_ == 128354),]
hist(ulig$rel_surfac)

junip <- data[which(data$taxref_cd_ == 104414),]
hist(junip$rel_surfac)


### plot % cover vs. altitude (priority > 1800m & > 80 %)
par(mfrow=c(2,2))

plot(rhodo$altinf, rhodo$rel_surfac, xlab="Elevation", ylab="% cover", main="Rhodo", pch=19,
     xlim=range(1000, 2600), ylim=range(0,100))
abline(v=1800, h=80, lty="dashed")

plot(myrt$altinf, myrt$rel_surfac, xlab="Elevation", ylab="% cover", main="Myrtille", pch=19,
     xlim=range(1000, 2600), ylim=range(0,100))
abline(v=1800, h=80, lty="dashed")

plot(ulig$altinf, ulig$rel_surfac, xlab="Elevation", ylab="% cover", main="Airelle", pch=19,
     xlim=range(1000, 2600), ylim=range(0,100))
abline(v=1800, h=80, lty="dashed")

plot(junip$altinf, junip$rel_surfac, xlab="Elevation", ylab="% cover", main="Juniper", pch=19,
     xlim=range(1000, 2600), ylim=range(0,100))
abline(v=1800, h=80, lty="dashed")


### bilan des courses : il nous faut + de relevés !! 


