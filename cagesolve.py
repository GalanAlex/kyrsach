import numpy as np
import matplotlib.pyplot as plt

import solve as slv

MODE = 'c'

J = 100
K = 100

Lx = 10
Ly = 6
L = 4
lambd = 2
n = 1
z = 1
x = 1

y = 3

o = complex(0,1)
k1 = 2*np.pi/lambd
alpha = o*(np.pi**2)/(4*k1*n*Lx**2)


balpha = np.zeros((J))
p=np.zeros((J,K))
l=np.zeros((J,K))





def cons():
    return((np.pi**2)/(Lx**2))
def psi(j,h):
    return np.exp(-(((j*h-Ly/2)/(0.1*Ly))**2))

def mainsolve(J,K):
    
    V = np.zeros((J, K))

    if MODE == 'c': 
        h = Ly/J #Hy
        tau = L/K #Hz
        betta = alpha/2

        a = abs(1 - tau*alpha - ((2*betta*tau)/(h**2)))
        b = abs(betta*tau/h**2)
        c = abs(1 + tau*alpha + ((2*betta*tau)/(h**2)))

        j=0
        while j< J-1:
            V[j,0]=psi(j,h)
            j=j+1
        print(V[:,0])

        balpha[0] = -b / a
        balpha[1] = b / (b * balpha[0] + a)
        m = 2
        while m <= J - 3:
            balpha[m] = -(b / (b * balpha[m - 1] + a))
            m = m + 1

        kk = 0
        while kk <=K-2:
            p[0, kk] = c * V[1, kk] - b * V[2, kk]
            p[J - 2, kk] = -b * V[J - 3, kk] + c * V[J - 2, kk]
            l[0, kk] = p[0, kk] * 1 / a
            m = 1
            while m <= J - 3:
                p[m, kk] = -b * V[m , kk] + c * V[m+1, kk] - b * V[m + 2, kk]
                l[m, kk] = (p[m, kk] - b * l[m - 1, kk]) / (b * balpha[m - 1] + a)
                V[J - 2, kk + 1] = (p[J - 2, kk] - b * l[J - 3, kk]) / b * balpha[J - 3] + a
                m = m + 1
            m = J-3
            while m != 1:
                V[m, kk + 1] = balpha[m] * V[m + 1, kk + 1] + l[m, kk]
                m = m - 1
            kk=kk+1
        #print(p[:,10])

        '''
        #print(a)
        #print(b)
        #print(c)
        
        j=0
        while j< J-1:
            V[j,0]=psi(j,h)
            print(V[j,0])
            j=j+1
        balpha[0] = -b / a
        balpha[1] = b / (b * balpha[0] + a)
        m = 2
        while m <= J - 3:
            balpha[m] = -(b / (b * balpha[m - 1] + a))
            m = m + 1
        
        #print(balpha)
        #print("-----------")
        p[0, 0] = c * V[1, 0] - b * V[2, 0]
        #print(p[0, 0])
        #print(V[1, 0])

        p[J - 2, 0] = -b * V[J - 3, 0] + c * V[J - 2, 0]
        #print(p[J - 2, 0])

        kk=0
        while kk <= K-2:
            m = 1

            while m <= J - 3:
                p[m, kk] = -b * V[m , kk] + c * V[m+1, kk] - b * V[m + 2, kk]
                #print(m," ",kk," ",p[m, kk])
                #print(V[m , kk])
                #print(V[m+1, kk])
                #print(V[m + 2, kk])

                m = m + 1
            kk=kk+1
        #print(p[:,3])
        #print()
        #print("------------")

        l[0, 0] = p[0, 0] * 1 / a
        kk=0
        while kk <= K-1:
            m = 1

            while m <= J - 3:
                l[m, kk] = (p[m, kk] - b * l[m - 1, kk]) / (b * balpha[m - 1] + a)
                m = m + 1
            kk=kk+1

        #print(l[25,25])
        #print("-----------")
        #print(V[25,0])
        V[J-2,kk+1]= (p[J-2,kk] - b*l[J-3,kk])/b*balpha[J-3]+a
        #print( V[J-3,1])
        #print("-----------")
        kk=0
        while kk !=K-2:
            m = J - 2

            while m != 1 :
                V[m,kk+1] = balpha[m]*V[m+1,kk+1]+l[m,kk]
                m=m-1
            kk=kk+1
        u=20
        while u <=20:
            #print("Интерация ",u,"для k  V - ",V[u,1]," p - ",p[u,1]," бета - ",l[u,1] )
            #print("Интерация ",u,"для k+1 V - ",V[u,2]," p - ",p[u,2]," бета - ",l[u,2])
            u=u+1
        '''
        '''
        A = np.zeros((J+1, J+1))
        B = np.zeros((J+1, J+1))
    
        j = 0 #заполнение матриц и начального условия
        while j < J-1:
            V[j, 0] = psi(j,h)
            A[j, j-1] = abs((betta*tau)/h**2)
            A[j, j] = abs(1 - tau*((alpha*np.pi**2)/(2*Lx**2)) - 2*betta*tau/h**2)
            A[j, j+1] = abs((betta*tau)/h**2)
            B[j, j-1] = abs((betta*tau)/h**2)
            B[j, j] = abs(1 + tau*((alpha*np.pi**2)/(2*Lx**2)) - ((2 * betta * tau )/( h**2)))
            B[j, j+1] = abs((betta*tau)/h**2)
            j = j+1

        V[0,0] = psi(0,h)
        V[J, 0] = psi(J, h)
        print(V)
        A[0,0] = abs(2*betta*tau/h**2)
        A[0,1] = abs(2*a*tau/(h**2))
        A[J, J] = 1 + a * tau / (h ** 2)
        A[J, J - 1] = -a * tau / (h ** 2)
        B[0,0] = abs(2*betta*tau/h**2)
        B[0,1] = abs(2*a*tau/(h**2))
        B[J, J] = 1 - a* tau / (h ** 2)
        B[J, J - 1] = a * tau / (h ** 2)

        k = 0
        while k < K:
            V[:, k+1] = np.linalg.solve(A, B.dot(V[:,k]))  #Решение для Av(k+1) = BV(k)
            k = k+1
        '''


    if MODE == 'o':
        h = Ly / J
        tau = L / K

        j = 0
        while j <= J-1:
            V[j, 0] = psi(j,h)
            j = j+1
        k = 0

        while k < K-1:
            V[0, k+1] = 0
            V[J-1, k + 1] = 0
            j = 1
            while j < J-1:
                V[j,k+1] = abs(alpha*tau*((np.pi**2)/(Lx**2))*V[j,k] - (alpha*tau*(V[j+1,k] - 2*V[j,k] + V[j-1, k]))/(h**2) + V[j,k])
                #print(V[j,k+1])
                j = j+1
            k = k+1



    if MODE == 'im':

        h = Ly / J  # Hy
        tau = L / K  # Hz
        #mu = alpha * tau/(h**2)
        #print(mu)
        beta = tau*alpha * ((np.pi**2)/(Lx**2))
        #print(beta)
        #betta = alpha / 2

        #A = np.zeros((J , J))

        mu = alpha * tau/(h**2)
        beta = tau*alpha * ((np.pi**2)/(Lx**2))



        #a = abs(1 - tau * alpha - ((2 * betta * tau) / (h ** 2)))
        #b = abs(betta * tau / h ** 2)
        #c = abs(1 + tau * alpha + ((2 * betta * tau) / (h ** 2)))
        a = abs(1-2*mu-beta)
        b = abs(mu)
        j = 0
        while j < J - 1:
            V[j, 0] = psi(j, h)
            #print(V[j,0])
            j = j + 1

        balpha[0] = -b / a
        balpha[1] = b / (b * balpha[0] + a)
        m = 2
        while m <= J - 3:
            balpha[m] = -(b / (b * balpha[m - 1] + a))
            m = m + 1

        kk = 0
        while kk <= K - 2:
            #p[0, kk] = c * V[1, kk] - b * V[2, kk]
            p[0,kk] = V[1,kk]
            #p[J - 2, kk] = -b * V[J - 3, kk] + c * V[J - 2, kk]
            p[J-2,kk]=V[J-1,kk]
            l[0, kk] = p[0, kk] / a
            m = 1
            while m <= J - 3:
                #p[m, kk] = -b * V[m, kk] + c * V[m + 1, kk] - b * V[m + 2, kk]
                p[m,kk] = V[m+1,kk]
                l[m, kk] = (p[m, kk] - b * l[m - 1, kk]) / (b * balpha[m - 1] + a)
                m = m + 1
            V[J - 1, kk + 1] = (p[J - 2, kk] - b * l[J - 3, kk]) / (b * balpha[J - 3] + a)

            m = J - 2
            while m != 1:
                V[m, kk + 1] = balpha[m-1] * V[m + 1, kk + 1] + l[m-1, kk]
                m = m - 1
            V[1,kk+1] = balpha[0]*V[2,kk+1]+l[0,kk]
            kk = kk + 1



        '''
        j = 0
        while j < J-1:
            V[j, 0] = psi(j,h)
            #print(psi(j,h))
            A[j+1, j] = abs(mu)
            #print(A[j+1, j])
            A[j, j] = abs(1 - 2 * mu - beta)
            #print( A[j, j])
            A[j, j+1] = abs(mu)
            j = j+1
        print(A)

        V[0,0] = psi(0,h)
        #print(psi(J,h))
        V[J-1, 0] = psi(J, h)
        #print(psi(J, h))
        A[0,0] = abs(1 - 2 * mu - beta)
        A[0,1] = abs(mu)
        A[J-1, J-1] = abs(1 - 2 * mu - beta)
        A[J-1, J-2 ] = abs(mu)
        #print(A.shape)
        


        k = 0
        while k < K-1:
            V[:, k+1] = np.linalg.solve(A, V[:,k])
            k = k+1
        '''
    return V



if __name__ == "__main__":
    (mainsolve(J,K))