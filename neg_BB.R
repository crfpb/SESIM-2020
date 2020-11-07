##############################################
#      food-web negative interactions        #
#  Firkowski et al. - Ecological Monographs  #
##############################################

# body sizes
P.aurelia    = (111.6) #non
P.bursaria   = (101.3) #non
Spirostomum  = (843.8) #non
Colpidium    = (81.0)  #prey
T.pyriformis = (26.7)  #prey
Vorticella   = (79.7)  #prey
Monostyla    = (127.3) #prey
Euplotes     = (85.4)  #pred
Blepharisma  = (471.3) #pred

# interactions
# effect of collumn species (.c) on row species (r.)

# constants
o = 0.0007
p = 0.75
e = 0.1

# P. aurelia
i4.4  = 0
i4.5  = 0
i4.6  = 0
i4.7  = 0
i4.8  = 0
i4.9  = 0
i4.10 = 0
i4.11 = 0
i4.12 = 0

# P. bursaria
i5.4  = 0
i5.5  = 0
i5.6  = 0
i5.7  = 0
i5.8  = 0
i5.9  = 0
i5.10 = 0
i5.11 = 0
i5.12 = 0

# Spirostomum sp.
i6.4  = 0
i6.5  = 0
i6.6  = 0
i6.7  = 0
i6.8  = 0
i6.9  = 0
i6.10 = 0
i6.11 = 0
i6.12 = 0

# Colpidium sp.
i7.4  = 0
i7.5  = 0
i7.6  = 0
i7.7  = 0
i7.8  = 0
i7.9  = 0
i7.10 = 0
i7.11 = -o*(Euplotes/Colpidium)^p
i7.12 = -o*(Blepharisma/Colpidium)^p

# Tetrahymena pyriformis
i8.4  = 0
i8.5  = 0
i8.6  = 0
i8.7  = 0
i8.8  = 0
i8.9  = 0
i8.10 = 0
i8.11 = -o*(Euplotes/T.pyriformis)^p
i8.12  = -o*(Blepharisma/T.pyriformis)^p

# Vorticella sp.
i9.4  = 0
i9.5  = 0
i9.6  = 0
i9.7  = 0
i9.8  = 0
i9.9  = 0
i9.10 = 0
i9.11 = -o*(Euplotes/Vorticella)^p
i9.12  = -o*(Blepharisma/Vorticella)^p

# Monostyla sp.
i10.4  = 0
i10.5  = 0
i10.6  = 0
i10.7  = 0
i10.8  = 0
i10.9  = 0
i10.10 = 0
i10.11 = 0
i10.12 = -o*(Blepharisma/Monostyla)^p

# Euplotes
i11.4  = 0
i11.5  = 0
i11.6  = 0
i11.7  = 0
i11.8  = 0
i11.9  = 0
i11.10 = 0
i11.11 = 0
i11.12 = -o*(Blepharisma/Euplotes)^p 

# Blepharisma
i12.4  = 0  
i12.5  = 0
i12.6  = 0
i12.7  = 0
i12.8  = 0
i12.9  = 0
i12.10 = 0
i12.11 = 0
i12.12 = -o*(Blepharisma/Blepharisma)^p

# food-web
fullFW_neg <- (rbind(cbind(i4.4,i4.5,i4.6,i4.7,i4.8,i4.9,i4.10,i4.11,i4.12),
                 cbind(i5.4,i5.5,i5.6,i5.7,i5.8,i5.9,i5.10,i5.11,i5.12),
                 cbind(i6.4,i6.5,i6.6,i6.7,i6.8,i6.9,i6.10,i6.11,i6.12),
                 cbind(i7.4,i7.5,i7.6,i7.7,i7.8,i7.9,i7.10,i7.11,i7.12),
                 cbind(i8.4,i8.5,i8.6,i8.7,i8.8,i8.9,i8.10,i8.11,i8.12),
                 cbind(i9.4,i9.5,i9.6,i9.7,i9.8,i9.9,i9.10,i9.11,i9.12),
                 cbind(i10.4,i10.5,i10.6,i10.7,i10.8,i10.9,i10.10,i10.11,i10.12),
                 cbind(i11.4,i11.5,i11.6,i11.7,i11.8,i11.9,i11.10,i11.11,i11.12),
                 cbind(i12.4,i12.5,i12.6,i12.7,i12.8,i12.9,i12.10,i12.11,i12.12)))
