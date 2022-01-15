def pbcDistAB(xA, yA, zA, xB, yB, zB, cellDim):
    from numpy import sqrt
     
    #impose PBCs
    #on x
    if abs(xA - xB) > cellDim["a"]/2:
        if xB < xA:   
            xB =   xB + cellDim["a"]
        else:
            xB =   xB - cellDim["a"]
    #on y
    if abs(yA - yB) > cellDim["b"]/2:
        if yB < yA:
            yB = yB + cellDim["b"]
        else:
            yB = yB - cellDim["b"]
    #on z
    #if abs(zA - zB) > cellDim["c"]/2:
    #    if zB < zA:
    #        zB = zB + cellDim["c"]
    #    else:
    #        zB = zB - cellDim["c"]
    
    xij = xA - xB
    yij = yA - yB
    zij = zA - zB
    rAB = sqrt(xij**2 + yij**2 + zij**2)
    
    return(xB, yB, zB, rAB)
