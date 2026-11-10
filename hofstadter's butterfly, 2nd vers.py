from matplotlib import pyplot as plt
import numpy as np
import seaborn as sns
import scipy as sc

plt.close('all')
a = 1
J = 2
i_c = complex(0,1)
alpha = 0.4

e_site = 1.5

def E(qx,qy): #band energy of the vector q
    alpha = -2*J
    c1 = np.cos(qx*a)
    c2 = np.cos(qy*a)
    Eq = alpha*(c1 + c2)

    return Eq

def norm2(V):
    S = 0
    n = len(V)
    for k in range(n):
        S += V[k]**2
    return np.sqrt(S)

def metric(v1,v2): #distance between 2 vectors, with the 2-norm
    V = v1-v2 #there are same lenght arrays
    metrique = norm2(V)
    return metrique


## 3d plot environment
x = np.linspace(-np.pi/a,np.pi/a,100)
y = np.linspace(-np.pi/a,np.pi/a,100)

X,Y = np.meshgrid(x,y)
"""
Z = E(X,Y)

fig = plt.figure()
ax = plt.axes(projection = '3d')
ax.set_facecolor('purple')
ax.plot_surface(X,Y,Z)
plt.show()
"""
##

#faut prendre une maille élémentaire tq, en raison du réseau carré, puisse être complètement divisé sans sortir de la maille (ie si on prend une maille de 3 par 1 ou l'inverse, faut prendre un nombre N = 3k, k entier naturel)


#en ce qui concerne la maille qu'on décide de conséidérer, on définira un array de la taille choisi avec les énergies d'intéraction entre les atomes bien choisi selon les atomes considérés dans la maille

"""
def produit_scal
def vect_plus_proche_voisin(R0,N,a1,a2):#on prend un vecteur d'une certaine taille et faut faire des tests vis-à-vis des coordonnées, s'il on sort de la maille ou pas, et ak sont les vecteurs élementaire de la maille en question
    i0,j0 = R0 #R0 est un élément de la matrice décrivant la maille 2D, de la forme (i0,j0)
    L_deltai =[]
    #faut trouver un code qui fonctionne pour chaque maille
    #regarder vis à vis des normes et prendre les plus petites qui seront égales entre elles

"""

## Class to define lattice (test)
class Polygone:
    def __init__(self,cote):
        self.cote=cote

    def nb_cote(self):
        return len(self.cote)

class Carre(Polygone):

    def __init__(self,c):
        super().__init__([c,c,c,c])
        self.c = c

    def vecteur_directeur(self):
        self.a1 = np.array([1,0])
        self.a2 = np.array([0,1])
        return self.a1, self.a2

class Trait:
    def __init__(self,point1,point2):
        self.point1 = point1
        self.point2 = point2

    def distance(self):
        self.d = metric(self.point1,self.point2)
        return self.d

## I'll see later for the hexagonal lattice

#définir, pour le dictionnaire contenant les points et les coordonnées, que l'on se place dans une maille de côté a (en unité de a) et donc que pour les coordonnées va falloir écrire *a, et pour les mailles va falloir renvoyer les vecteurs qui forment le réseau de Bravais

### LATTICE LENGHT

R0 = np.array([0,0])
Nx = 10
Ny = 10
def lattice_dico(Nx,Ny): #the function will only build the dictionnary for the squared lattice
    lattice = {}

    lx = np.arange(0,Nx,1,dtype = np.int64)  #makes the research of dots easier than putting some float numbers
    ly = np.arange(0,Ny,1,dtype = np.int64)

    X,Y = np.meshgrid(lx,ly)

    for i in range(len(lx)):
        for j in range(len(ly)):
            Xij = X[i,j]
            Yij = Y[i,j]
            Rij = (Xij,Yij) #I really need to obtain all differents dots as a key to get the whole lattice bc of dictionaries constraints
            Rij_calc_norm = np.array([Xij,Yij]) #see above
            norm_Rij = norm2(Rij_calc_norm)
            lattice[Rij]=norm_Rij
    return lattice

## CAUTION ! Board effects

def nearest_neighbours(R,lattice): #R is an array (to avoid all doubts)
    x,y = R

    Rkey = (x,y)
    lattice_copy = lattice
    normR = lattice.pop(Rkey) #we gain the norm and it avoid a 0 in the list of distances

    norm_difference = {}
    L_nearest_neighbours = []
    L_norm_diff = []
    for dot in lattice.keys():
        xd,yd = dot
        Rd = np.array([xd,yd])
        dist_from_R = metric(Rd,R)

        norm_difference[(xd,yd)] = dist_from_R #to simplify the research

        L_norm_diff.append(dist_from_R)

    #for the squared lattice, we know the result and if we draw a circle around the R point with a radius equal to the minimal distance of a nearest neigbours, we can see that all nearest neighbours are on the sphere (mathematical definition of "sphere"). I don't know if this results works on others lattices but i guess i can use it here.
    rmin = min(L_norm_diff)

    for dott in norm_difference.keys():
        Rij = dott
        distij = norm_difference[Rij]

        if distij <= rmin:
            L_nearest_neighbours.append(Rij)

    return L_nearest_neighbours

"""
R1 = np.array([0,7])

lattice_sq = lattice_dico(Nx,Ny)

L_nn = nearest_neighbours(R1,lattice_sq)

print(L_nn)

"""

def delta_nearest_n(R,L_nn):

    L_delta = []
    x,y = R
    n = len(L_nn)
    for k in range(n):
        Rk = L_nn[k]
        xk,yk = Rk

        x_delta = x - xk

        y_delta = y - yk

        deltas = np.array([x_delta,y_delta])

        L_delta.append(deltas)
    return L_delta
"""
R1 = np.array([0,7])

lattice_sq = lattice_dico(Nx,Ny)

L_nn = nearest_neighbours(R1,lattice_sq)

print(L_nn)

deltas = delta_nearest_n(R1,L_nn)
print(deltas)
"""
### Energy spectrum

# To compute interesting vallues, we have to take N>>1: by using tight biding approximation, Wave functions are Bloch functions with several orbitals
# So, The Hofstadter's butterfly will depend of the number of bands that we consider, and also the number of orbitals per elementary lattices.

#To start: create the matrix of interactions.

#board conditions (didnt noticed that before)

def Hamiltonian(N,J,E_site,alpha): #Only matrices in Mn(C) can be diagonalized
    Areal = np.zeros((N,N))
    Aimag = np.zeros_like(Areal)
    epsilon_J = E_site/(-J) #because we would multiply the whole matrix by tunneling coefficient J at the end of the process
    diagonal_X = [np.ones(N),epsilon_J*np.ones(N),np.ones(N)]
    offsets = np.array([-1,0,1])
    Hamilt_magnetic_less = sc.sparse.dia_matrix((diagonal_X,offsets),shape=(N,N))
    L_a = [Areal,Aimag]
    Hamiltonian_magnetic_less = Hamilt_magnetic_less.toarray()
    for A in L_a:
        for j in range(N):
            for l in range(N):
                A[j,l] = Hamiltonian_magnetic_less[j,l]

        for j in range(1,N):
            A[0,j] = 1
            A[N-1,j-1] = 1
            A[j,0] = 1
            A[j-1,N-1] = 1 # to respect boundary conditions

        initial_phase = i_c*2*np.pi*alpha


    for k in range(2):
        for l in range(N):
            for j in range(N):
                Ak = L_a[k]
                phase_l = initial_phase*j
                phase_term = np.exp(phase_l)

                Ak[j,l] = Ak[j,l]*((1-k)*phase_term.real + k*phase_term.imag)
    #compute in the "for" loop real and imaginary part, python explodes and solutions are NOT EASY to understand

    return -J*Areal, -J*Aimag

## HAMILTONIAN CREATION

Hr, Hi = Hamiltonian(Nx,J,e_site,alpha)

#print(Hr,Hi)

H = Hr + i_c*Hi


Hcg = np.conjugate(H)
Hcg_ = np.zeros_like(Hcg)

for i in range(N):
    for j in range(N):
        Hcg_[i,j] = Hcg[j,i]

H_final = H + Hcg_ #to get an hermitian hamiltonian and to be able to diagonalize it

print(H_final)


## H_final is the total hamiltonian of the system, hermitian...

##############################################################################

## creation of the hamiltonian for a elementary lattice:

def solutions_elementary_lattice(L_R_tuple,p_alpha,q):#here, as alpha = p'/p (to simplify the analysis), and in calculations show that only p appears in the results (p and p' are integers). L_R_tuple contains the coordinates of the x points of the elementary lattice, q is the Bloch momentum that would be useful to compute bands.
    return None





























