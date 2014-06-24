from __future__ import division
from math import *
import numpy as np

def ikJacobian2DOF(Length1, Length2, alpha1, alpha2, oriX, oriY,  x, y):
    maxIterations=300
    L1=Length1
    L2=Length2
    x=float(x)
    y=float(y)
    

    theta1=alpha1/(180/pi)
    theta2=alpha2/(180/pi)
    
    x0=oriX+L1*cos(theta1)+L2*cos(theta1-theta2)
    y0=oriY+L1*sin(theta1)+L2*sin(theta1-theta2)

    iDist=sqrt((x-x0)**2+(y-y0)**2)

    #print("Current Distance:"+str(round(iDist,2))+"xRot:"+str(int(xRot))+"yRot:"+str(int(yRot)))
    x1=(x-x0)/1
    y1=(y-y0)/1

    if iDist>0:

        Jacobian = [[0,0],[0,0]]
        Jacobian[0][0]=-L1*sin(theta1)-L2*sin(theta1-theta2)
        Jacobian[0][1]=L2*sin(theta1-theta2)
        Jacobian[1][0]=L1*cos(theta1)+L2*cos(theta1-theta2)
        Jacobian[1][1]=-L2*cos(theta1-theta2)

        det=Jacobian[0][0]*Jacobian[1][1]-Jacobian[0][1]*Jacobian[1][0]
        if (det==0):
            det=1.0
        
        invJacobian = [[0,0],[0,0]]
        invJacobian[0][0]=Jacobian[1][1]/det
        invJacobian[0][1]=-Jacobian[0][1]/det
        invJacobian[1][0]=-Jacobian[1][0]/det
        invJacobian[1][1]=Jacobian[0][0]/det

        deltaAlpha1=invJacobian[0][0]*(x1)+invJacobian[0][1]*(y1)
        deltaAlpha2=invJacobian[1][0]*(x1)+invJacobian[1][1]*(y1)

        alpha1+=deltaAlpha1
        alpha2+=deltaAlpha2
         
    return (deltaAlpha1, deltaAlpha2)

def ikAnalytical3DOF(coxaLength, femurLength, tibiaLength, oriXYZ, destXYZ):
    ##
    ## Based on Oscar Liang tutorial : http://blog.oscarliang.net/inverse-kinematics-and-trigonometry-basics/
    ## Takes a XYZ position and calculates the angles need to position angular joints of hexapod limb
    ## Note: Angles are different from Oscar's. The relation:
    ## My Angels  |  Oscar's Angles
    ##   alpha    >    gamma         (Coxa Angle)
    ##   beta     >    alpha - 90      (FemurAngle in relation to Coxa)
    ##   gamma    >    180 - beta      (TibiaAngle in relation to Femur)
    
    #print('Moving to:'+str(destXYZ))
    xf=float(destXYZ[0])
    yf=float(destXYZ[1])
    Zoffset=float(oriXYZ[2])-float(destXYZ[2])


    L1=sqrt(xf**2+yf**2)
    #print('L1='+str(L1))

##  Changed gamma angle computation from original since I work with angles (+90/-90) around X axis (not Y)
    if xf==0:
        if yf>0:
            gamma=pi/2
        else:
            gamma=-pi/2
    else:
        gamma=atan(yf/xf) 

   
    L=sqrt(Zoffset**2+(L1-coxaLength)**2)

    try:
        alpha1=acos(Zoffset/L)   
        alpha2=acos((tibiaLength**2-femurLength**2-L**2)/(-2*femurLength*L)) ##multiplicar racio por -1 para simplificar
        alpha=alpha1+alpha2

        beta=acos((L**2-tibiaLength**2-femurLength**2)/(-2*tibiaLength*femurLength)) ## "    
    except Exception, e:
        print('Error: Unable to calculate angles...', e)
        return False
    else:
        pass
    finally:
        pass

    ## Proceeding to conversion acoording to description
    alphaNew=gamma
    betaNew =alpha-pi/2
    gammaNew=pi-beta

    #print('alpha='+str(degrees(alphaNew)))
    #print('beta='+str(degrees(betaNew)))
    #print('gamma='+str(degrees(gammaNew)))
        
    return (alphaNew, betaNew, gammaNew)

def ikJacobian3DOF(L1, L2, L3, alpha, beta, gamma, origin, targetPosition, maxAngleDisp, precision=0.5, maxiterations=100):
    #alpha, beta and gamma in radians.
    listofangularinc=[] # Returns a list of all necessary angular displacements

    def calcPos(L1, L2, L3, alpha, beta, gamma, origin):
        position=[0,0,0]
        
        projL=L1+L2*cos(beta)+L3*cos(beta-gamma)
        position[0]=origin[0] + cos(alpha)*projL
        position[1]=origin[1] + sin(alpha)*projL
        position[2]=origin[2] + L2*sin(beta) + L3*sin(beta-gamma)

        return position
    
    def calcDist(targetPos, actualPos):
##        print('Distance to target: %s' % startDistance)
        return sqrt((targetPosition[0]-actualPos[0])**2+(targetPosition[1]-actualPos[1])**2+(targetPosition[2]-actualPos[2])**2)

    ## Starting Position:
    startposition=calcPos(L1, L2, L3, alpha, beta, gamma, origin)

    currentposition=startposition
    currentdistance=calcDist(targetPosition, currentposition)

    i=0
    while currentdistance > precision and i < maxiterations:
        i += 1

        jacobian=[[0,0,0],
                  [0,0,0],
                  [0,0,0]]

    ##    x = cos(alpha)*(L1+L2cos(beta)+L3cos(beta-gamma))    
    ##    y = sin(alpha)*(L1+L2cos(beta)+L3cos(beta-gamma)) 
    ##    z = L2sin(beta)+L3sin(beta-gamma)
        
                                                                         
        jacobian[0][0] = -sin(alpha)*(L1+L2*cos(beta)+L3*cos(beta-gamma))
        jacobian[0][1] = -cos(alpha)*(L2*sin(beta)+L3*sin(beta-gamma))
        jacobian[0][2] = L3*cos(alpha)*sin(beta-gamma)
        
        jacobian[1][0] = cos(alpha)*(L1+L2*cos(beta)+L3*cos(beta-gamma))
        jacobian[1][1] = -sin(alpha)*(L2*sin(beta)+L3*sin(beta-gamma))
        jacobian[1][2] = L3*sin(alpha)*sin(beta-gamma)
        
        jacobian[2][0] = 0
        jacobian[2][1] = L2*cos(beta)+L3*cos(beta-gamma)
        jacobian[2][2] = -L3*cos(beta-gamma)                           

        varPos=np.matrix([[targetPosition[0]-currentposition[0]],
                          [targetPosition[1]-currentposition[1]],
                          [targetPosition[2]-currentposition[2]]])
        
        ##varPos=varPos.T ##Confirmar!!!
        
        jacobianM=np.matrix(jacobian)

        pseudoInverseM=np.linalg.pinv(jacobianM)
        varAngles=pseudoInverseM*varPos

        currMaxAngle = abs(varAngles).max() # Returns the maximum absolute angle from the angles matrix

        if currMaxAngle > maxAngleDisp:
            fact = float(maxAngleDisp) / currMaxAngle

            varAngles = varAngles*fact

        (alpha, beta, gamma)=(alpha+varAngles[0,0], beta+varAngles[1,0], gamma+varAngles[2,0])

        # Updating current distance to target
        currentposition=calcPos(L1, L2, L3, alpha, beta, gamma, origin)
        currentdistance=calcDist(targetPosition, currentposition)

        listofangularinc.append((varAngles[0,0], varAngles[1,0], varAngles[2,0])) # Adds calculated angular displacements

        jacobian=None
        varPos=None
        jacobianM=None
        pseudoInverseM=None
        varAngles=None

    ##    print("varPos:",varPos)
    ##    print("Jacobian:",jacobian)
    ##    print("alpha:",degrees(varAngles[0]),"beta:",degrees(varAngles[1]),"gamma:",degrees(varAngles[2]))

    ##    print('ReachedPosition: ' ,calcPos(L1, L2, L3, varAngles[0], varAngles[1], varAngles[2], origin))
        ##print('JacobianInverse\n',jacobian.I)

    if currentdistance <= precision:
        return listofangularinc
    else:
        print('Unable to reach target position. Current distance:', currentdistance)
        return False

def getjointanglesforpath(coxalength, femurlength, tibialength, origin, p0, p1, nbrsteps=50):
    # Returns a list of tuples containing the sequencial joint angles (in degrees) to move the
    # limb along a linear path
    listofjointangles=[]
    p0=np.array(p0)
    p1=np.array(p1)
    convertradstointdegrees = lambda r: round(degrees(r),1)

    for i in range(1, nbrsteps+1):     
        t=i/nbrsteps
        #print(t)
        pm=linearbezier(p0, p1, t)
        #print('Moving to ', pm)
        jointangles=ikAnalytical3DOF(coxalength, femurlength, tibialength, origin, list(pm))
        listofjointangles.append(map(convertradstointdegrees, jointangles))

    return listofjointangles


def linearbezier(p0, p1, t):
    ## p0 and p1 are of type numpy.array
    return p0+t*(p1-p0)
