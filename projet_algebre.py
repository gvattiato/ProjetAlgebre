import numpy as numpy
import matplotlib as plt
import pylab as pl
# run a file on python:
# python projet.py


# Importation of the 2, 3 columns from the BIOGRID file
data = numpy.genfromtxt('BIOGRID-ORGANISM-Rattus_norvegicus-3.2.120.tab2.txt',skip_header=1,usecols = (1, 2))
print "data[:,0]=",data[:,0]
print "data[:,1]=",data[:,1]

#---------------------------------------------------------------
# Matrice d'adjacence A
#---------------------------------------------------------------
# We first need to identify how many proteins there are
# Combine the two columns of data
allgenes = numpy.concatenate([data[:,0],data[:,1]])
# commande: unique ????
# remove the duplicates and order the vector
genes = list(set(allgenes))
genes.sort()
#print "genes=",genes
print "Number of Genes ={0}".format(len(genes))
N=0
# Create a square matrix full of 0
A = numpy.zeros((len(genes),len(genes)))
# Fill up A
for i in range(0,len(data)):
	for j in  range(0,len(genes)):
		if(data[i,0] == genes[j]):
			for k in range(0,len(genes)):
				if(data[i,1] == genes[k]):
					A[j,k]+=1
					#print "A[",j,",",k,"]=",A[j,k]
					N+=1
print "Interactions: {0} and len(data)={1}".format(N,len(data))
print "A ={0}".format(A)


#------------------------------------------
A = numpy.genfromtxt('matrice_A.txt')




#---------------------------------------------------------------
# Mesures de centralite (score)
#---------------------------------------------------------------
# Vector with the degree of each protein
cd = numpy.zeros((len(genes),1))
# Index of the most central protein Cd
index_cd_max=0
for i in range(0,len(genes)):
	for j in range(0,len(genes)):
		cd[i]+=A[i,j]/(len(genes)-1)
	if( cd[i]>cd[index_cd_max] ):
		index_cd_max = i

print "degre de centralite:{0}".format(cd)
print "La proteine la plus centrale est {0}, de degre Cd={1}".format(genes[index_cd_max],cd[index_cd_max])



#---------------------------------------------------------------
# Calcul de la centralite par valeurs propres
#---------------------------------------------------------------



# Calcul des valeurs propres de A
eig_val_A, eig_vec_A = numpy.linalg.eig(A)
print "linalg.eig(A):\nvecteurs propres {v}\nvaleurs propres {w}".format(v = eig_vec_A, w = eig_val_A)
print "(On a bien {0} valeurs propres)".format(len(eig_val_A))
index_ce_max = 0 # Index of the most central protein Ce 
# Computation of the index of the highest eigen value of A
Ce = eig_vec_A[0]
print "Ce = {0}".format(Ce)
for i in range(0,len(Ce)):
	if (Ce[i] == 0):
		if( (Ce[i].real)>(Ce[index_ce_max].real)):
#		if(numpy.linalg.norm(eig_vec_A[i])>numpy.linalg.norm(eig_vec_A[index_ce_max])):
			index_ce_max = i
			print index_ce_max
print "La proteine la plus centrale est {0}, de Ce = {1}".format(genes[index_ce_max],Ce[index_ce_max])


#LSQ = numpy.linalg.lstsq(cd,(eig_vec_A).real)
#pl.plot(LSQ)

#---------------------------------------------------------------
# Comparaison de Cd et Ce par une régression linéaire
#---------------------------------------------------------------





#---------------------------------------------------------------
# Calcul des coefficients a et b et deuxième régression linéaire
#---------------------------------------------------------------





#---------------------------------------------------------------
# Ranking
#---------------------------------------------------------------







