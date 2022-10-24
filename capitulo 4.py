import libreria as lc
import vectores as lc
def amplitude_ket(n1, n2, k1, k2):
    k3 = k2
    r1 = lc.producto_Interno(n1, k1)
    r2 = lc.producto_Interno(n2, k2)
    return lc.producto_Interno(r1, r2)


def proba(h, ket):
    r = lc.norma(ket)
    x = ket[h]
    y = x
    sum_1 = lc.cplxmult(x, y)
    sum_2 = abs(sum_1[0] + sum_1[1])
    pro = sum_2/(r ** 2)
    return pro



def varianza(M, v):
    if len(M) == len(v[0]):
        H = lc.multmat(M, lc.trama(v))
        w = [[]]
        for j in H:
            w[0].append(j[0])
        x = lc.multmat(lc.trama(w), v)
        E = x[0][0][0] + x[1][0][1]
        m1 = lc.ident(E, M)
        N = lc.cplxrest(M, m1)
        Delta = lc.multmat(N, N)
        r = v
        for i in range(len(v)):
            for j in range(len(v[0])):
                x = v[i][j]
                c = x[0] ** 2
                t = x[1] ** 2
                r[i][j] = (c, t)
        Rf = lc.multmat(r, Delta)
        x = lc.prm(Rf[0][0])
        return round(x[0], 1)


def estadospropios(ob, v):
    ob1 = []
    for i in range(len(ob)):
        ob1.append([])
        for j in range(len(ob[0])):
            ob1[i].append(ob[i][j])
    for i in range(len(ob1)):
        for j in range(len(ob1[0])):
            ob1[i][j] = lc.complejo(ob1[i][j])
    N = list(ob1[0])
    y = varianza(ob, v)
    return N, y

def detectarParticula(ket,X):
    normaKet=lc.normaMatriz([ket])
    normaComplejo=lc.normaMatriz([[ket[X]]])
    return (normaComplejo**2)/(normaKet**2)
def transitarVector(ket,ket2):
    ket1=lc.transpuesta([ket])
    bra=lc.matrizConjugada([ket2])
    car=lc.multiplicacionMatrizMatriz(bra,ket1)[0][0]
    normaBra=lc.normaMatriz(bra)
    normaket=lc.normaMatriz(ket1)
    car2=(normaBra*normaket,0)
    return lc.division(car,car2)

def interfaz():
    ket=[tuple(map(float, x.split(","))) for x in (input("ket ").split(" "))]
    opt = int(input("Opcion"))
    while opt!=0:
        if (opt==1):
            print("Ingrese posicion")
            X= int(input("X: "))
            print(detectarParticula(ket,X))
        elif (opt==2):
            print("Ingrese vector ket")
            ket2=[tuple(map(float, x.split(","))) for x in (input("ket ").split(" "))]
            print(transitarVector(ket,ket2))
        else:
            print("Opcion Invalida")
        opt = int(input("Opcion "))

def ParticulaPosicion(ket,X):

    return lc.detectarParticula(ket,X)
def transitarVectorVector(ket,ket2):

    return lc.transitarVector(ket,ket2)
def valorEsperado(obs,ket):
    obsSobreket = lc.accion(obs, ket)
    bra = lc.matrizConjugada(obsSobreket)
    ket1 = lc.transpuesta([ket])
    bra1 = lc.transpuesta(bra)
    car = lc.multiplicacionMatrizMatriz(bra1, ket1)[0][0]
    return car

def varianzaObservable(obs,ket):
    nve = lc.multiplicacion(valorEsperado(obs, ket), (-1, 0))
    mve = lc.multiplicacionEscalarMatriz(lc.matrizIdentidad(len(obs)), nve)
    delta = lc.sumaMatrices(obs, mve)
    deltaCuadrado = lc.multiplicacionMatrizMatriz(delta, delta)
    var = valorEsperado(deltaCuadrado, ket)
    return var
def propiosObservable(obs):
    for i in range(len(obs)):
        for j in range(len(obs[0])):
            obs[i][j]=complex(obs[i][j][0],obs[i][j][1])
    return obs

def probabilidadObservable(obs,ket):

    valP,vectP = propiosObservable(obs)
    probs=[]
    for v in vectP:
        p=lc.transitarVector(v,ket)
        probs.append(p)
    return probs
def dinamica(un,init,steps):
    up=[un]
    ur=un
    for p in range(steps):
        ur=lc.multiplicacionMatrizMatriz(ur,ur)
        up.append(ur)
    un1=[]
    for k in range(len(up)):
        un1=lc.multiplicacionMatrizMatriz(up[k],up[k-1])
    temp=lc.accion(un1,init)
    return temp