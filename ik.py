from __future__ import division
import math
import numpy as np

def ikJacobian2DOF(Length1, Length2, alpha1, alpha2, oriX, oriY,  x, y):
    maxIterations=300
    L1=Length1
    L2=Length2
    x=float(x)
    y=float(y)
    

    theta1=alpha1/(180/math.pi)
    theta2=alpha2/(180/math.pi)
    
    x0=oriX+L1*math.cos(theta1)+L2*math.cos(theta1-theta2)
    y0=oriY+L1*math.sin(theta1)+L2*math.sin(theta1-theta2)

    iDist=math.sqrt((x-x0)**2+(y-y0)**2)

    #print("Current Distance:"+str(round(iDist,2))+"xRot:"+str(int(xRot))+"yRot:"+str(int(yRot)))
    x1=(x-x0)/1
    y1=(y-y0)/1

    if iDist>0:

        Jacobian = [[0,0],[0,0]]
        Jacobian[0][0]=-L1*math.sin(theta1)-L2*math.sin(theta1-theta2)
        Jacobian[0][1]=L2*math.sin(theta1-theta2)
        Jacobian[1][0]=L1*math.cos(theta1)+L2*math.cos(theta1-theta2)
        Jacobian[1][1]=-L2*math.cos(theta1-theta2)

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
    
    print('Moving to:'+str(destXYZ))
    xf=float(destXYZ[0])
    yf=float(destXYZ[1])
    Zoffset=float(oriXYZ[2])-float(destXYZ[2])


    L1=math.sqrt(xf**2+yf**2)
    print('L1='+str(L1))

##  Changed gamma angle computation from original since I work with angles (+90/-90) around X axis (not Y)
    if xf==0:
        if yf>0:
            gamma=math.pi/2
        else:
            gamma=-math.pi/2
    else:
        gamma=math.atan(yf/xf) 

   
    L=math.sqrt(Zoffset**2+(L1-coxaLength)**2)
    print('L='+str(L))

    alpha1=math.acos(Zoffset/L)   
    alpha2=math.acos((tibiaLength**2-femurLength**2-L**2)/(-2*femurLength*L)) ##multiplicar racio por -1 para simplificar
    alpha=alpha1+alpha2

    beta=math.acos((L**2-tibiaLength**2-femurLength**2)/(-2*tibiaLength*femurLength)) ## "


    ## Proceeding to conversion acoording to description
    alphaNew=gamma
    betaNew =alpha-math.pi/2
    gammaNew=math.pi-beta

    print('alpha='+str(math.degrees(alphaNew)))
    print('beta='+str(math.degrees(betaNew)))
    print('gamma='+str(math.degrees(gammaNew)))
        
    return (alphaNew, betaNew, gammaNew)

def ikJacobian3DOF(L1, L2, L3, alpha, beta, gamma, origin, targetPosition):
    #alpha, beta and gamma in radians.
    def calcPos(L1, L2, L3, alpha, beta, gamma, origin):
        position=[0,0,0]
        projL=L1+L2*math.cos(beta)+L3*math.cos(beta-gamma)
        position[0]=origin[0] + math.cos(alpha)*projL
        position[1]=origin[1] + math.sin(alpha)*projL
        position[2]=origin[2] + L2*math.sin(beta) + L3*math.sin(beta-gamma)

        return position

    
    def calcDist(targetPos, actualPos):
        startDistance=math.sqrt((targetPosition[0]-actualPos[0])**2+(targetPosition[1]-actualPos[1])**2+(targetPosition[2]-actualPos[2])**2)
##        print('Distance to target: %s' % startDistance)
        return startDistance

    ## Starting Position:
    startPosition=calcPos(L1, L2, L3, alpha, beta, gamma, origin)

    startDistance=calcDist(targetPosition, startPosition)

    jacobian=[[0,0,0],
              [0,0,0],
              [0,0,0]]

##    x = cos(alpha)*(L1+L2cos(beta)+L3cos(beta-gamma))    
##    y = sin(alpha)*(L1+L2cos(beta)+L3cos(beta-gamma)) 
##    z = L2sin(beta)+L3sin(beta-gamma)
    
                                                                     
    jacobian[0][0] = -math.sin(alpha)*(L1+L2*math.cos(beta)+L3*math.cos(beta-gamma))
    jacobian[0][1] = -math.cos(alpha)*(L2*math.sin(beta)+L3*math.sin(beta-gamma))
    jacobian[0][2] = L3*math.cos(alpha)*math.sin(beta-gamma)
    
    jacobian[1][0] = math.cos(alpha)*(L1+L2*math.cos(beta)+L3*math.cos(beta-gamma))
    jacobian[1][1] = -math.sin(alpha)*(L2*math.sin(beta)+L3*math.sin(beta-gamma))
    jacobian[1][2] = L3*math.sin(alpha)*math.sin(beta-gamma)
    
    jacobian[2][0] = 0
    jacobian[2][1] = L2*math.cos(beta)+L3*math.cos(beta-gamma)
    jacobian[2][2] = -L3*math.cos(beta-gamma)                           

    varPos=np.matrix([[targetPosition[0]-startPosition[0]],
                      [targetPosition[1]-startPosition[1]],
                      [targetPosition[2]-startPosition[2]]])
    
    ##varPos=varPos.T ##Confirmar!!!
    
    jacobianM=np.matrix(jacobian)
##    print(jacobian)
    
##    varAngles=jacobian*jacobian.T+LAMBDA**2*matrixI
##    varAngles=jacobian.T*varAngles.I
##    varAngles=varAngles*varPos

    pseudoInverseM=np.linalg.pinv(jacobianM)
    varAngles=pseudoInverseM*varPos

    
##    print("varPos:",varPos)
##    print("Jacobian:",jacobian)
##    print("alpha:",math.degrees(varAngles[0]),"beta:",math.degrees(varAngles[1]),"gamma:",math.degrees(varAngles[2]))

##    print('ReachedPosition: ' ,calcPos(L1, L2, L3, varAngles[0], varAngles[1], varAngles[2], origin))
    ##print('JacobianInverse\n',jacobian.I)
         
    return varAngles
