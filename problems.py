from scipy import special
from fdscheme import *

class Problem1:
    def __init__(self, I = 1.0):
        self.I = I # the inertia number

        self.yinf = 5
        self.ny   = 51
        self.y  = np.linspace(0, self.yinf, self.ny)
        self.dy = self.y[1] - self.y[0]

        self.tend = 5
        self.nt = self.tend*200+1
        self.t = np.linspace(0, self.tend, self.nt)
        self.dt = self.t[1] - self.t[0]

        print('Time and space resolution: ', self.dt, self.dy)
        
        self.us = 1
        self.ts= np.zeros_like(self.t)
        self.ts[0] = self.us # time series of us
        self.u = np.zeros_like(self.y)
        self.u[0] = 1
        
    def solve(self):
        u = self.u
        un = np.zeros_like(self.y)
        
        for n in range(1, self.nt):
            un = u + laplacianCar(self.y, u)*self.dt

            # update boundary condition
            un[0] = self.us
            un[-1]= un[-2]
            
            # update solid velocity
            dudy0 = (-1.5*u[0]+2*u[1]-0.5*u[2])/self.dy
            self.us += 1/self.I*dudy0*self.dt
            self.ts[n] = self.us

            u[:] = un[:]

class Problem2:
    def __init__(self, I = 1.0):
        self.I = I # the inertia number

        self.yinf = 1
        self.ny   = 21
        self.y  = np.linspace(0, self.yinf, self.ny)
        self.dy = self.y[1] - self.y[0]

        self.tend = 5
        self.nt = self.tend*1000+1
        self.t = np.linspace(0, self.tend, self.nt)
        self.dt = self.t[1] - self.t[0]

        print('Time and space resolution: ', self.dt, self.dy)
        
        self.us = 1
        self.ts= np.zeros_like(self.t)
        self.ts[0] = self.us # ts of us
        self.u = np.zeros_like(self.y)
        self.u[0] = 1
        
    def solve(self):
        u = self.u
        un = np.zeros_like(self.y)
        
        for n in range(1, self.nt):
            un = u + laplacianCar(self.y, u)*self.dt
            
            # update solid velocity
            dudy0 = (-1.5*u[0]+2*u[1]-0.5*u[2])/self.dy
            self.us += 1/self.I*dudy0*self.dt
            self.ts[n] = self.us

            # update boundary condition
            un[0] = self.us
            un[-1]= 0

            # update time
            u[:] = un[:]

class Problem3:
    def __init__(self, I = 1.0):
        self.I = I

        self.rinf = 10
        self.nr = 101
        self.r = np.linspace(1, self.rinf, self.nr)
        self.dr = self.r[1] - self.r[0]

        self.tend = 12
        self.nt = self.tend*1000+1
        self.t = np.linspace(0, self.tend, self.nt)
        self.dt = self.t[1] - self.t[0]

        print('Time and space resolution: ', self.dt, self.dr)
        self.us = 1
        self.ts= np.zeros_like(self.t)
        self.ts[0] = self.us # ts of us
        self.u = np.zeros_like(self.r)
        self.u[0] = 1

    def solve(self):
        u = self.u
        un = np.zeros_like(self.r)

        for n in range(1, self.nt):
            # semi-implicit update
            rhs = laplacianCyl(self.r, u)
            un = (u/self.dt + rhs)/(1/self.dt + 1/(self.r**2))

            # update solid velocity
            dudr0 = (-3*u[0]+4*u[1]-u[2])/(2*self.dr) - u[0]
            self.us += (4.0/self.I*dudr0)*self.dt
            self.ts[n] = self.us

            # update boundary condition
            un[0] = self.us
            un[-1] = un[-2]
            
            # update time
            u[:] = un[:]

class Problem4:
    def __init__(self, I = 1.0):
        self.I = I

        self.rinf = 10
        self.nr = 101
        self.r = np.linspace(1, self.rinf, self.nr)
        self.dr = self.r[1] - self.r[0]

        self.tend = 12
        self.nt = self.tend*1000+1
        self.t = np.linspace(0, self.tend, self.nt)
        self.dt = self.t[1] - self.t[0]

        print('Time and space resolution: ', self.dt, self.dr)
        self.us = 1
        self.ts= np.zeros_like(self.t)
        self.ts[0] = self.us # ts of us
        self.u = np.zeros_like(self.r)
        self.u[0] = 1

    def solve(self):
        u = self.u
        un = np.zeros_like(self.r)

        for n in range(1, self.nt):
            # semi-implicit update
            rhs = laplacianSph(self.r, u)
            un = (u/self.dt + rhs)/(1/self.dt + 2/(self.r**2))

            # update solid velocity
            dudr0 = (-3*u[0]+4*u[1]-u[2])/(2*self.dr) - u[0]
            self.us += (5.0/self.I*dudr0)*self.dt
            self.ts[n] = self.us

            # update boundary condition
            un[0] = self.us
            un[-1] = un[-2]
            
            # update time
            u[:] = un[:]

def CheckP2():
    print('Checking problem 2')
    pa = Problem2(0.25)
    pb = Problem2(1)
    pc = Problem2(4)
    pa.solve()
    pb.solve()
    pc.solve()

    mma_data = np.loadtxt("mma_p2.dat")

    plt.figure(figsize = (3.50394*1.4,1.4*2.16555))
    plt.semilogy(pa.t, pa.ts,'r-',label='I=1/4')
    plt.semilogy(pb.t, pb.ts,'g-',label='I=1'  )
    plt.semilogy(pc.t, pc.ts,'b-',label='I=4'  )

    plt.semilogy(mma_data[:,2], mma_data[:,3], 'go')
    plt.semilogy(mma_data[:,0], mma_data[:,1], 'ro')
    plt.semilogy(mma_data[:,4], mma_data[:,5], 'bo')

    # the universal function from P1
    plt.semilogy(pa.t, np.exp(pa.t)*special.erfc(np.sqrt(pa.t)), 'k:')

    rootsqa = 0.230491
    rootsqb = 0.740174
    rootsqc = 1.59919
    plt.semilogy(pa.t,0.9208*np.exp(-pa.t*rootsqa),'k--')
    plt.semilogy(pb.t,0.7299*np.exp(-pb.t*rootsqb),'k--')
    plt.semilogy(pc.t,0.3704*np.exp(-pc.t*rootsqc),'k--')

    plt.tight_layout()
    plt.xlim(0, 2)
    plt.ylim(1e-2, 1.1)
    plt.xlabel('t')
    plt.ylabel(r'$u_s$')
    plt.savefig('p2.svg')
    plt.show()

def CheckP3P4():
    print('Checking problem 3 and 4')
    pa = Problem3(0.25)
    pb = Problem3(1)
    pc = Problem3(4)
    pa.solve()
    pb.solve()
    pc.solve()

    myt = np.linspace(1,12)
    myu = 1/32.0*myt**(-2)

    mma_data = np.loadtxt("mma_p34.dat")
    plt.figure(figsize = (3.50394*1.4,1.4*2.16555))
    plt.loglog(pa.t, pa.ts*4,'r-',label='I=1/4')
    plt.loglog(pb.t, pb.ts,'g-',label='I=1'  )
    plt.loglog(pc.t, pc.ts/4,'b-',label='I=4'  )
    plt.loglog(mma_data[:,6], 4*mma_data[:,7], 'ro')
    plt.loglog(mma_data[:,8], mma_data[:,9], 'go')
    plt.loglog(mma_data[:,10],0.25*mma_data[:,11],'bo')
    plt.loglog(myt, myu, 'k-', lw=2)

    print('Checking problem 4')
    pa = Problem4(0.25)
    pb = Problem4(1)
    pc = Problem4(4)
    pa.solve()
    pb.solve()
    pc.solve()
    myu = 1/(60*np.sqrt(np.pi))*myt**(-2.5)
    plt.loglog(pa.t, pa.ts*4,'r--',label='I=1/4')
    plt.loglog(pb.t, pb.ts,'g--',label='I=1'  )
    plt.loglog(pc.t, pc.ts/4,'b--',label='I=4'  )
    plt.loglog(mma_data[:,12], mma_data[:,13]*4, 'r*')
    plt.loglog(mma_data[:,14], mma_data[:,15], 'g*')
    plt.loglog(mma_data[:,16], mma_data[:,17]/4, 'b*')
    plt.loglog(myt, myu, 'k-', lw=2)

    plt.tight_layout()
    plt.xlim(0.1, 12)
    plt.ylim(1e-5, 1)
    plt.xlabel('$t$')
    plt.ylabel(r'$u_s/I$')
    plt.savefig('p3_4.svg')
    plt.show()

if __name__ == "__main__":
    CheckP2()
    CheckP3P4()
