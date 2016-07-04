'''modelling system of a trebuchet
assumptions
- rigid frame
- constant torque friction at the pivot proportional to wieght
- no air resistance while before release
- release happens before weight stall point
- main weight before release hangs straight down
- mass of main arm is uniformly distributed
- slings are always under tension
- relative to the weight and arm the rotation of the sling is insignificant

Algebra;

State of trebuchet variables;
a = angle of main arm
b = angle of firing sling
da = change in a (angular velocity)
db = change in b (angular velocity)
Ts = tension in sling
Ta = torque turning arm
Fb = force turning sling
KE = kinetic energy 
GPE = gravitational potential energy
FE = frictional energy
T = tension in projectile sling
x = actual distance of projectile to pivot
y = actual distance of weight to pivot

Design of trebuchet;
m = small mass (projectile)
M = large mass (firing weight)
RM, Rm, Rl = radius from pivot along arm
rM, rm = radius of swing from end of arm
l = lenght of arm
d = density of arm
h = hieght of pivot above slider
'''

#init conditions: design
#D = dict(RM = 1, Rm = 4, Rl = (Rm - RM)/2, rM = 1, rm = 4, M = 100, m = 1, d = 1, g = 9.81, h = 3)

from tkinter import *
from math import *
from time import *

def build (dt0=0.01,a=1):
    global l, RM, Rm, Rl, rM, rm, M, m, d, g, h, j0, dt, F
    RM = 1
    Rm = 4
    Rl = (Rm - RM)/2
    l = Rm + RM
    rM = 1
    rm = 4.5
    M = 25
    m = 1
    d = 1
    g = 9.81
    h = 4
    j0 = (a,acos((cos(a)*Rm - h)/rm),0,0,0,0,0)
    dt = dt0
    #0.01 is the relative frictional coefficient
    F = 0.01*(M+m+l*d)
    

#energy conservation section
def GPE (a,b):
    return (cos(a) * (RM*M -l*d*Rl - m*Rm) + cos(a+b)*rm*m)*g

def Inertia (a,b,c=0):
    # about the pivot
    global x, y
    x = Rm - cos(b)*rm
    y = RM - cos(a+c)*rM
    return m*x**2 + M*y**2 + l*d*Rl**2

def calcda (gpe,I,fe,initgpe,db):
    return sqrt(2*(initgpe - fe - gpe)/I)

#newtonian calculations for main swinging arm
def calcdda (a,da,b,db,I):
    return (sin(a)*M*g*RM - l*d*Rl*sin(a) - sin(b)*Rm*(db**2*rm*m - cos(a+b)*m*g-cos(b)*da**2*Rm*m)-F)/(I + (sin(b)*Rm)**2*m)

def calcddaslide(a,da,b,I):
    return (sin(a)*M*g*RM - l*d*Rl*sin(a))/(I + (sin(b)*Rm)**2*m)
    
#newtonian calculations of projectile sling
def calcddb (a,da,dda,b):
    return (g*sin(a+b) + da**2*Rm*sin(b)  + cos(b)*dda*Rm)/rm

#calculates a single time jump
def timejump (a,b,da,db,dda,ddb,t):
    global slide
    if slide:
        I1 = Inertia(a,b)
        a += da*dt
        b = acos((cos(a)*Rm - h)/rm)-a
        I2 = Inertia(a,b)
        da = da*I1/I2 + dda*dt
        dda = calcddaslide(a,da,b,I2)
        t+=dt
        if dda*sin(b)*cos(a+b)*Rm > g or a+b >= pi:
            slide = False
            db = (sin(a)*da*Rm)/(sin(b)*rm)
            print('off slide')
    else:
        I1 = Inertia(a,b)
        a += da*dt
        b += db*dt
        I2 = Inertia(a,b)
        da = da*I1/I2 + dda*dt
        db += ddb*dt
        dda = calcdda(a,da,b,db,I2)
        ddb = calcddb(a,da,dda,b)
        t+=dt
    return a,b,da,db,dda,ddb,t

def release(a,b,da,db):
    ang1 = 2*pi -(a+b)
    x = Rm - cos(b)*rm
    ang2 = pi - a + asin(rm*sin(b)/x)
    v = [(db*rm*cos(ang1)+da*x*cos(ang2)),(db*rm*sin(ang1)+da*x*sin(ang2))]   
    return v

#intiates loop that fires the trebuchet
def fire ():
    global slide, angledata
    # slide is whether the shot it on the slide
    slide = True
    last = 10/dt
    n=0
    j=j0
    angledata=[]
    #initgpe = GPE(j[0],j[1])
    while n < last and j[1] < 2.7:
        n+=1
        if (n-1)%(int(last/200)) == 0:
            angledata += [[round(ans,3) for ans in j]]
        j = timejump(j[0],j[1],j[2],j[3],j[4],j[5],j[6])
    angledata += [[round(ans,3) for ans in j]]
    print(angledata)
    print(release(j[0],j[1],j[2],j[3]))
    #runge kutta method needed for later version

def coordadd(i,j):
    return [round(i[0] + j[0]),round(i[1] + j[1])]

# from angles it calculates the coordinates of parts of the trebuchet
def coordcalc (a,b):
    #needs to be scaled to size
    scale = 20
    centerpivot = [200,350-scale*h]
    slingpivot = coordadd(centerpivot,[-Rm*sin(a)*scale,Rm*cos(a)*scale])
    weightpivot = coordadd(centerpivot,[RM*sin(a)*scale,-RM*cos(a)*scale])
    shot = coordadd(slingpivot,[rm*sin(a+b)*scale,-rm*cos(a+b)*scale])
    weight = coordadd(weightpivot,[0,scale*rM])
    return [centerpivot, slingpivot, weightpivot, shot, weight]

build()   
fire()

coords = coordcalc(0.5,acos((cos(0.5)*Rm - h)/rm))

def animate(n=0):
    global coords
    #coordsdata = [coordcalc(j[0],j[1]) for j in angledata]
    j = angledata[n]
    coords = coordcalc(j[0],j[1])
    print(coords)
    #for x in coords:
        #C.create_circle(x[0],x[1],5)
    if n+1 < len(angledata):
        #sleep(0.5)
        for x in coords:
            C.create_oval(x[0]-10,x[1]-10,x[0]+10,x[1]+10,fill='red')
        animate(n+1)
        



top = Tk()
L1 = Label(top, text="rm")
L1.pack(side = TOP)
E1 = Entry(top, bd=5, textvariable = StringVar())
E1.pack()
B2 = Button(top, text="animate", command=animate)
B2.pack()
B1 = Button(top,text="fire", command=fire)
B1.pack()
C = Canvas(top, width = 400, height = 400)
C.pack(side = RIGHT)
for x in coords:
    C.create_oval(x[0]-10,x[1]-10,x[0]+10,x[1]+10,fill='red')

top.mainloop()












