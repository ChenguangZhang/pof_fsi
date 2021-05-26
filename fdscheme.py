import numpy as np
import matplotlib.pyplot as plt

def deltau(u):
    # central difference
    return np.r_[-3*u[0]+4*u[1]-u[2],  u[2:]-u[:-2],  3*u[-1]-4*u[-2]+u[-3]]

def laplacianCyl(r, u):
    # (1/r) (r u_r)_r
    dr = r[1] - r[0]
    dudr = deltau(u)/(2*dr)
    return deltau(r*dudr)/(2*dr)/r

def laplacianSph(r, u):
    # (1/r*r) (r*r u_r)_r
    dr = r[1] - r[0]
    dudr = deltau(u)/(2*dr)
    return deltau(r*r*dudr)/(2*dr)/(r*r)

def laplacianCar(y, u):
    # u_yy
    dy = y[1] - y[0]
    return np.r_[2*u[0]-5*u[1]+4*u[2]-u[3], u[2:]+u[:-2]-2*u[1:-1], 2*u[-1]-5*u[-2]+4*u[-3]-u[-4]]/(dy*dy)

def test_sine():
    x = np.linspace(0, 2*np.pi)
    y = np.sin(x)
    plt.plot(x, laplacianCar(x,y), 'r--o', x, -np.sin(x), 'k-')
    plt.show()

if __name__ == "__main__":
    test_sine()
