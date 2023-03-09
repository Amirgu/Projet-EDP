from typing import List
import gmsh
import numpy as np
import random as rdm
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.path import Path
import matplotlib.patches as patches
import matplotlib.tri as mtri


def LoadVTX(s):
    x = []
    y = []
    i = 0
    with open(s, "r") as filin:
        for l in filin:
            coord = l.split()
            if 5 == len(coord):
                x.append([int(coord[1]), int(coord[2]), int(coord[3])])
                y.append(int(coord[4]))

                i = i + 1
    x = np.array(x)
    y = np.array(y)
    return x, y


def LoadELT(s):
    a = []
    b = []
    with open(s, "r") as filin:
        for l in filin:
            co = l.split()
            if len(co) == 3:
                a.append(float(co[1]))
                b.append(float(co[2]))
    return a, b


x = LoadVTX("config1.txt")[0] ### ELEMENTS
a = LoadELT("config1.txt")[0] ### VTX
b = LoadELT("config1.txt")[1] ### VTX
y = LoadVTX("config1.txt")[1] ### Label
a = np.asarray(a)
b = np.asarray(b)


def plotdomains(s):                                                 ### Q2 : en fonction du nom du fichier
                                                                    ### je l'ai fait pour faciliter l'appel
    x = LoadVTX(s)[0]
    y = LoadVTX(s)[1]
    a = LoadELT(s)[0]
    b = LoadELT(s)[1]
    color = ["blue", "red"]
    for i in range(len(x)):
        if y[i] == 1:

            u = [a[x[i][0]], a[x[i][1]], a[x[i][2]], a[x[i][0]]]
            v = [b[x[i][0]], b[x[i][1]], b[x[i][2]], b[x[i][0]]]
            plt.plot(u, v, color[0])
            plt.fill(u, v, color='blue', alpha=.76)
            plt.title("OMEGA1 ET OMEGA 2")
        else:
            u = [a[x[i][0]], a[x[i][1]], a[x[i][2]]]
            v = [b[x[i][0]], b[x[i][1]], b[x[i][2]]]
            plt.plot(u, v, color[1])
            plt.fill(u, v, color='red', alpha=0.76)

    plt.axis("equal")
    plt.xlim(-2, 2)


    red_patch = patches.Patch(color='red',label='Omega2')
    blue_patch = patches.Patch(color='blue',label='Omega1')
    plt.legend(handles=[red_patch,blue_patch])

    plt.savefig("proxy_artists.png")

    plt.show()




def prodscal(v, u):
    s = 0
    for i in range(len(v)):
        s += v[i] * u[i]

    return s


def airetri(j):                                     ###calcul l'aire d'un triangle j donnée dans le tableau des vtx
    s1 = [a[x[j][0]], b[x[j][0]]]
    s2 = [a[x[j][1]], b[x[j][1]]]
    s3 = [a[x[j][2]], b[x[j][2]]]
    l1 = [s1[0] - s2[0], s1[1] - s2[1]]
    l2 = [s2[0] - s3[0], s2[1] - s3[1]]
    l3 = [s3[0] - s1[0], s3[1] - s1[1]]
    T = np.array([[s2[0] - s1[0], s2[1] - s1[1]], [s3[0] - s1[0], s3[1] - s1[1]]])
    u = np.linalg.det(T)
    u = 0.5 * abs(u)

    return u


def airT(u, v, w):
    s1 = [u[0], u[1]]
    s2 = [v[0], v[1]]
    s3 = [w[0], w[1]]
    l1 = [s1[0] - s2[0], s1[1] - s2[1]]
    l2 = [s2[0] - s3[0], s2[1] - s3[1]]
    l3 = [s3[0] - s1[0], s3[1] - s1[1]]
    T = np.array([[s2[0] - s1[0], s2[1] - s1[1]], [s3[0] - s1[0], s3[1] - s1[1]]])
    u = np.linalg.det(T)
    u = 0.5 * abs(u)
    return u


def moins(x, y):                            #### Fonction utile
    s = np.zeros(2)
    s[0], s[1] = x[0] - y[0], x[1] - y[1]
    return s


def prodv(s):                               #### Fonction utile 2
                                            #### utilisé dans le calul du vecteur normal d'un triangle
    v = np.zeros(2)
    v[0], v[1] = s[1], -s[0]
    return v


def multp(l, v):
    v = np.multiply(v, l)
    return v


def val(j, k, l):                               #### Calcul la rigidité élementaire pour un triangle l
                                                #### en fonction de j,k dans [0,2]
    s1 = np.array([a[x[l][0]], b[x[l][0]]])
    s2 = np.array([a[x[l][1]], b[x[l][1]]])
    s3 = np.array([a[x[l][2]], b[x[l][2]]])
    f = np.array([s1, s2, s3, s1, s2, s3])
    s = prodscal(moins(f[j + 2], f[j + 1]), moins(f[k + 2], f[k + 1]))
    s = s / (4 * airetri(l))
    return s


def lagloc(l, k, t1, t2):
    s1 = np.array([a[x[l][0]], b[x[l][0]]])
    o = 0                                                   ##### fonction de forme locale pour un noeud
    s2 = np.array([a[x[l][1]], b[x[l][1]]])                 ##### pour un noeud i si il est dans le triangle l
    s3 = np.array([a[x[l][2]], b[x[l][2]]])                 ##### et qu'il est le k-eme element du triangle
    f = np.array([s1, s2, s3, s1, s2, s3])                  ##### t1,t2 sont les points calculés
    v1 = prodv(s1)
    v2 = prodv(s2)
    v3 = prodv(s3)
    n = np.array([v1, v2, v3])
    verts = [
        (a[x[l][0]], b[x[l][0]]),  # left, bottom
        (a[x[l][1]], b[x[l][1]]),  # left, top
        (a[x[l][2]], b[x[l][2]]),  # right, top
        (0., 0.),  # right, bottom
        # ignored
    ]

    codes = [
        Path.MOVETO,
        Path.LINETO,
        Path.LINETO,
        Path.CLOSEPOLY,
    ]

    path = Path(verts, codes)
    jrt = (t1, t2)
    flag = path.contains_point(jrt)                                     #### utilisée pour savoir si un point t1,t2 est dans un triangle

    if flag == True:
        o = (prodscal(moins(np.array([t1, t2]), f[k]), n[k]))
        o = o /(prodscal(moins(f[k + 1], f[k]), n[k]))
        o = 1 - o
    return o


def lagfunct(i, t1, t2):
    s = np.array((len(x), 1))
    v = 0
    for j in range(len(x)):
        k = 0
        while k < 3:
            verts = [
                (a[x[j][0]], b[x[j][0]]),  # left, bottom
                (a[x[j][1]], b[x[j][1]]),  # left, top
                (a[x[j][2]], b[x[j][2]]),  # right, top
                (77., 0.),  # right, bottom
                # ignored
            ]
            codes = [
                Path.MOVETO,
                Path.LINETO,
                Path.LINETO,
                Path.CLOSEPOLY,
            ]

            path = Path(verts, codes)
            jrt = (t1, t2)
            flag = path.contains_point(jrt)


            if x[j][k] == i :
                ax = plt.plot()
                patch = patches.PathPatch(path, facecolor='orange', lw=2)
                ax.add_patch(patch)
                ax.set_xlim(-2, 2)
                ax.set_ylim(-2, 2)

                print(lagloc(j,k,t1,t2))


            k += 1

    plt.show()


def val2(j, k, l):                                 ### Pour calculer la contribution d'un élement dans une triangle
                                                   ### au lieu de calculer la matrice de la contribution d'un élement
    if y[l] == 2:                                  ### pour la matrice de masse et Omega 2
        if j == k:
            return airetri(l) / 6
        else:
            return airetri(l) / 12
    else:
        return 0

def val3(j,k,l):

    if j == k:
        return airetri(l) / 6
    else:
        return airetri(l) / 12


def estdanslebord(j):
    u = x[j]                                            #### test un triangle si il est dans le bord du domaine
    t1 = u[0]                                           #### et renvoie 1(vrai,faux),les position (,) pour savoir ou
    t2 = u[1]                                           #### renvoie aussi la longeur de l'arete qui est dans le bord
    t3 = u[2]
    if (b[t1] == 0.25 and b[t2] == 0.25):
        l=moins([a[t1], b[t1]], [a[t2], b[t2]])
        longueur=np.sqrt(prodscal(l,l ))
        return 1,(0,1),2,longueur
    if (b[t2] == 0.25 and b[t3] == 0.25):
        l = moins([a[t3], b[t3]], [a[t2], b[t2]])
        longueur = np.sqrt(prodscal(l, l))
        return 1,(1,2),2,longueur
    if (b[t1] == 0.25 and b[t3] == 0.25):
        l = moins([a[t1], b[t1]], [a[t3], b[t3]])
        longueur = np.sqrt(prodscal(l, l))
        return 1,(0,2),2,longueur
    if (b[t1] == -0.25 and b[t2] == -0.25):
        l = moins([a[t1], b[t1]], [a[t2], b[t2]])
        longueur = np.sqrt(prodscal(l, l))
        return 1,(0,1),2,longueur
    if (b[t2] == -0.25 and b[t3] == -0.25):
        l = moins([a[t3], b[t3]], [a[t2], b[t2]])
        longueur = np.sqrt(prodscal(l, l))
        return 1, (1, 2),2,longueur
    if (b[t1] == -0.25 and b[t3] == -0.25):
        l = moins([a[t1], b[t1]], [a[t3], b[t3]])
        longueur = np.sqrt(prodscal(l, l))
        return 1,(0,2),2,longueur
    if (a[t1] == -0.25 and a[t2] == -0.25):
        l = moins([a[t1], b[t1]], [a[t2], b[t2]])
        longueur = np.sqrt(prodscal(l, l))
        return 1, (0, 1),1,longueur
    if (a[t2] == -0.25 and a[t3] == -0.25):
        l = moins([a[t3], b[t3]], [a[t2], b[t2]])
        longueur = np.sqrt(prodscal(l, l))
        return 1, (1, 2),1,longueur
    if (a[t1] == -0.25 and a[t3] == -0.25):
        l = moins([a[t1], b[t1]], [a[t3], b[t3]])
        longueur = np.sqrt(prodscal(l, l))
        return 1, (0, 2),1,longueur
    if (a[t1] == -0.75 and a[t2] == -0.75):
        l = moins([a[t1], b[t1]], [a[t2], b[t2]])
        longueur = np.sqrt(prodscal(l, l))
        return 1, (0, 1),1,longueur
    if (a[t2] == -0.75 and a[t3] == -0.75):
        l = moins([a[t3], b[t3]], [a[t2], b[t2]])
        longueur = np.sqrt(prodscal(l, l))
        return 1, (1, 2),1,longueur
    if (a[t1] == -0.75 and a[t3] == -0.75):
        l = moins([a[t1], b[t1]], [a[t3], b[t3]])
        longueur = np.sqrt(prodscal(l, l))
        return 1, (0, 2),1,longueur

    return 0,0


def longueurbord(A, u):                   #####Pour tester l'exactitude de la methode
                                          ##### Renvoie normalement 2
    n = len(A)
    s = 0
    v = []
    for i in range(n):
        s = 0
        for j in range(n):
            s += A[i][j] * u[j]

        v.append(s)
    return v


def finalv():
    y = len(b)
    A = np.zeros((y, y))
    n = len(x)
    s = 0
    for i in range(n):
        for m in range(3):
            for u in range(3):
                j = x[i][m]
                k = x[i][u]
                s = val(m, u, i) + val2(m, u, i)            #### ajout des contributions élementaires
                A[j][k] += s
                if estdanslebord(i)[0]:                     #### si dans le bord
                    if (u == m and (estdanslebord(i)[1][0] == m or estdanslebord(i)[1][1] == m)):           ### rajoute les coefficients de masse de bord
                            A[j][k]+= estdanslebord(i)[3] / 3                                               ### je justifierais de la formule dans un annexe
                    elif  (u!= m) and ((estdanslebord(i)[1][0] == m and estdanslebord(i)[1][1] == u) or (estdanslebord(i)[1][0] == u and estdanslebord(i)[1][1] == m)):
                            A[j][k]+= estdanslebord(i)[3] / 6                   ####estdanslebord[3] longueur segment



    return A


def finalB():
    B = np.zeros((len(b), 1))

    for i in range(len(x)):
        for k in range(3):
            if y[i] == 2:                                           ### si dans omega 2
                j = x[i][k]
                B[j] += (airetri(j)) / 3                             ### rajoute le coefficient de second membre
            if estdanslebord(i)[0]:                                             #### si dans le bord
                if estdanslebord(i)[1][0] == k:
                    l=estdanslebord(i)[1][1]
                    t1=a[x[i][k]]+a[x[i][l]]                                    ### calcul de la quadrature directement
                    t1=t1/2
                    t2 = b[x[i][k]] + b[x[i][l]]
                    t2 = t2 / 2
                    B[j] += ((estdanslebord(i)[3]) / 6) * (1 + lagloc(i, k, t1, t2))
                if estdanslebord(i)[1][1] == k:
                    l=estdanslebord(i)[1][0]
                    t1=a[x[i][k]]+a[x[i][l]]
                    t1=t1/2
                    t2 = b[x[i][k]] + b[x[i][l]]
                    t2 = t2 / 2
                    B[j] += ((estdanslebord(i)[3]) / 6) * (1 + lagloc(i, l, t1, t2))





    return B

def rechertr():
    s = []
    v = []
    u = []
    p=[]
    t=0
    for i in range(len(x)):
        if estdanslebord(i)[0]:

            s.append(i)

            u = [ a[ x[i][0] ], a[ x[i][1] ], a[ x[i][2] ], a[ x[i][0] ] ]

            v = [ b[ x[i][0] ], b[ x[i][1] ], b[ x[i][2] ], b[ x[i][0] ] ]

            p.append((estdanslebord(i)[3], i))
            t+=estdanslebord(i)[3]
            plt.plot(u,v)



A=finalv()
B=finalB()
X=np.linalg.solve(A,B)


def erreurquadratique():        #####Q6:calcul direct de l'érreur
    s=0
    for i in range(len(x)):
        for m in range(3):
            for u in range(3):
                j=x[i][m]
                k=x[i][u]
                s+=val3(m,u,i)*(X[j]-1)*(X[k]-1)        ### val3: c'est pour calculer l'integrale des formes en un triangle dans tout OMEGA

    s=np.sqrt(s)
    return s                                        ### erreur quadratique à peu prés 1/2 ce qui est enorme
                                                    ### ce la peut s'explique par la non-linearité de l'équation


def finalv2():
    y = len(b)
    A1 = np.zeros((y, y))
    n = len(x)
    s = 0
    for i in range(n):
        for m in range(3):
            for u in range(3):
                j = x[i][m]
                k = x[i][u]
                s = val(m, u, i) + val3(m, u, i)  #### Ajout de val3 à la place de val2 car c'est ce qui change dans la formulation variationelle
                A1[j][k] += s
                if estdanslebord(i)[0]:  #### si dans le bord
                    if (u == m and (estdanslebord(i)[1][0] == m or estdanslebord(i)[1][
                        1] == m)):  ### rajoute les coefficients de masse de bord
                        A1[j][k] += estdanslebord(i)[3] / 3  ### je justifierais de la formule dans un annexe
                    elif (u != m) and ((estdanslebord(i)[1][0] == m and estdanslebord(i)[1][1] == u) or (
                            estdanslebord(i)[1][0] == u and estdanslebord(i)[1][1] == m)):
                        A1[j][k] += estdanslebord(i)[3] / 6

    return A1

def finalB2():
    B1 = np.zeros((len(b), 1))

    for i in range(len(x)):
        for k in range(3):
            if y[i] == 2:                                           ### si dans omega 2
                j = x[i][k]
                B1[j] += (airetri(j)) / 3                             ### rajoute le coefficient de second membre                                                        ### supression du terme de bord dans le second membre
    return B1

def moyetvar(X):
    s = 0                                               ####erreur moyenne et ecart type
    for i in range(len(X)):
        s += X[i]-1

    s = s / len(X)
    t = 0
    for i in range(len(X)):
        t += (s - X[i]) * (s - X[i])
    t = t / len(X)
    t=np.sqrt(t)
    return(s,t)                                                     #### moyenne à 0.95 et ecart type à 0.4
def Graphe2D(S):
    triangles = x
    triang = mtri.Triangulation(a, b, triangles)

    z = np.zeros((len(a)))
    for i in range(len(z)):

            z[i]=S[i]


    t = plt.tricontourf(triang, z)

    plt.title('Graphe de la solution en 2D',
              fontsize=14, fontweight='bold')
    plt.axis("equal")
    plt.xlim(-1 - 0.05, 1 + 0.05)
    plt.ylim(-0.5 - 0.005, 0.5 + 0.005)

    plt.colorbar(t)

    plt.show()

def Graphe3D(S):
    triangles = x
    triang = mtri.Triangulation(a, b, triangles)
    T1 = np.linspace(-1,1,len(a))
    T2 = np.linspace(-0.5,0.5,len(b))

    for i in range(len(T1)):

            T1[i],T2[i]=a[i],b[i]


    z = np.sin(-T1 * T2)
    for i in range(len(z)):

            z[i]=S[i]-1


    ax = plt.figure().add_subplot(projection='3d')

    ax.plot_trisurf(T1, T2,triangles ,z, linewidth=0.3, antialiased=True)
    plt.title('Solution en 3D de la deuxieme equation')
    plt.show()
A1=finalv2()
B1=finalB2()
X1=np.linalg.solve(A1,B1)
Graphe3D(X)
def testmasse():
    y = len(b)
    Y = np.zeros((y, y))
    n = len(x)
    s = 0
    for i in range(n):
        for m in range(3):
            for u in range(3):
                j = x[i][m]
                k = x[i][u]
                s = val(m, u, i)
                Y[j][k] += s

    U=np.ones(len(b))
    Y=np.dot(Y,U)
    print(Y)
    t=0
    for i in range(len(Y)):
        t+=Y[i]
    print(t/len(Y))                                         #### Nous donne une valeur en e-17 ce qui est normale avec les erreurs des arrondis


def testrigidité():
    y = len(b)
    Y = np.zeros((y, y))
    n = len(x)
    s = 0
    for i in range(n):
        for m in range(3):
            for u in range(3):
                j = x[i][m]
                k = x[i][u]
                s = val3(m, u, i)
                Y[j][k] += s

    U=np.ones(len(b))
    T=np.dot(Y,U)
    print(prodscal(T,U))               #### Nous donne 1.75 ce qui est attendu pour val3 la surface de tout omega
                                       #### pour avoir la surface de Omega 2 il suffit de change val3 par val2




