'''modelling system of a trebuchet
assumptions
- rigid frame
- constant torque friction at the pivot proportional to wieght
- no air resistance
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

from tkinter import *
from math import *
from time import *

def build (dt0=0.01,a=0.5):
    global l, RM, Rm, Rl, rM, rm, M, m, d, g, h, j0, dt, F, releaseangle
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
    # corrects error if sling doesn't react the slide
    if (cos(a)*Rm - h)/rm > -1:
        j0 = (a,acos((cos(a)*Rm - h)/rm)-a,0,0,0,0,0)
    else:
        j0 = (a,pi-a,0,0,0,0,0)
    dt = dt0
    #0.01 is the relative frictional coefficient
    F = 0.01*(M+m+l*d)
    # sets the release angle for the sling
    releaseangle = 2.7

def Inertia (a,b,c=0):
    # calculates  inertia about the pivot
    global x, y
    x = Rm - cos(b)*rm
    y = RM - cos(a+c)*rM
    return m*x**2 + M*y**2 + l*d*(l**2/12+Rl**2) # correct main arm

#newtonian calculations for main swinging arm to get the angular acceleration
def calcdda (a,da,b,db,I):
    return (sin(a)*M*g*RM - l*d*Rl*sin(a) - sin(b)*Rm*(db**2*rm*m - cos(a+b)*m*g-cos(b)*da**2*Rm*m)-F)/(I + (sin(b)*Rm)**2*m)

def calcddaslide(a,da,b,I):
    #several terms removed as physics while the shot is on the slide
    return (sin(a)*M*g*RM - l*d*Rl*sin(a))/(I + (sin(b)*Rm)**2*m)
    
#newtonian calculations of projectile sling to get the angular acceleration
def calcddb (a,da,dda,b):
    return (g*sin(a+b) + da**2*Rm*sin(b)  + cos(b)*dda*Rm)/rm

#calculates a single time jump
def timejump (a,b,da,db,dda,ddb,t):
    global slide
    # this does the time jump for when the shot is on the slide
    if slide:
        I1 = Inertia(a,b)
        a += da*dt
        # conditions for end of slide
        # for a short while one of these conditions was a+b >= 3
        if dda*sin(b)*cos(a+b)*Rm > g  or ((cos(a)*Rm - h)/rm) < - 1:
            slide = False
            db = (sin(a)*da*Rm)/(sin(b)*rm)
        else:
            b = acos((cos(a)*Rm - h)/rm)-a
        I2 = Inertia(a,b)
        #conservation of angular moment
        da = da*I1/I2 + dda*dt
        dda = calcddaslide(a,da,b,I2)
        t+=dt
        
    else:
        #and then when it moves off the slide
        I1 = Inertia(a,b)
        a += da*dt
        b += db*dt
        I2 = Inertia(a,b)
        #conservation of angular moment
        da = da*I1/I2 + dda*dt
        db += ddb*dt
        dda = calcdda(a,da,b,db,I2)
        ddb = calcddb(a,da,dda,b)
        t+=dt
    return a,b,da,db,dda,ddb,t

#calculates the release velocity of the shot
def release(a,b,da,db):
    ang1 = 2*pi -(a+b)
    x = Rm - cos(b)*rm
    ang2 = pi - a + asin(rm*sin(b)/x)
    v = [(db*rm*cos(ang1)+da*x*cos(ang2)),(db*rm*sin(ang1)+da*x*sin(ang2))]   
    return v

#intiates loop that fires the trebuchet
def fire ():
    global slide, angledata, v
    # slide is whether the shot it on the slide
    slide = True
    last = 10/dt
    n=0
    j=j0
    angledata=[]
    while n < last and j[1] < releaseangle:
        n+=1
        #time conditions to record the data ever 0.05 seconds
        if (n-1)%(int(last/200)) == 0:
            angledata += [[round(ans,3) for ans in j]]
        j = timejump(j[0],j[1],j[2],j[3],j[4],j[5],j[6])
    angledata += [[round(ans,3) for ans in j]]
    #print(angledata)
    v = (release(j[0],j[1],j[2],j[3]))
    #runge kutta method needed for later version

# simply adds two vectors
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
    return [scale,centerpivot, slingpivot, weightpivot, shot, weight]

# runs the current design
build()   
fire()

coords = coordcalc(0.5,acos((cos(0.5)*Rm - h)/rm))

# animates the firing of the trebuchet
def animate(n=0):
    global coords
    while n+1 < len(angledata):
        sleep(0.1)
        draw(n)
        n+=1

# draw each part of the trebuchet
def draw(n):
    j = angledata[n]
    x = coordcalc(j[0],j[1])
    #0=scale,1=centre,2=slingpivot,3=weightpivot,4=shot,5=weight
    #print(coords)
    C.delete('all')
    # frame
    C.create_line(200-x[0]*2,350,x[1][0],x[1][1],fill='brown',width=7)
    C.create_line(200+x[0]*2,350,x[1][0],x[1][1],fill='brown',width=7)
    #main arm
    C.create_line(x[2][0],x[2][1],x[3][0],x[3][1],fill='black',width=5)
    #counterweight rod
    C.create_line(x[5][0],x[5][1],x[3][0],x[3][1],fill='black',width=4)
    #sling
    C.create_line(x[4][0],x[4][1],x[2][0],x[2][1],fill='black',width=2)
    #counterweight
    C.create_rectangle(x[5][0]-10,x[5][1]-10,x[5][0]+10,x[5][1]+10,fill='grey')
    C.update()


#small function used to tune the trebuchet
def drange(start, stop, step):
    r = start
    while r < stop:
        yield r
        r += step

#prints the data of the angles, angular velocity, angular acceleration and time
def printangles():
    print(angledata)

#prints velocity
def printvelocity(): 
    if v[1] > 0 and v[0] > 0:
        statement = 'Good release'
    elif v[1] > 0:
        statement = 'Early release'
    else:
        statement = 'Late release'
    print(v,statement)

#creates the gui
top = Tk()
B2 = Button(top, text="animate", command=animate)
B2.pack()
B1 = Button(top,text="print angle data", command=printangles)
B1.pack()
B3 = Button(top,text="print release velocity", command=printvelocity)
B3.pack()
C = Canvas(top, width = 400, height = 400)
C.pack(side = RIGHT)

top.mainloop()










