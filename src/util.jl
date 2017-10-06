
function GAMMAH(N::Integer, over::Real)
    #INTEGER N,I,J
    #DOUBLE PRECISION OVER
    I = N
    J = 2 * I - 1
    retval = 1.0
    while (J>=1) && retval < over
         retval = retval * J/2.0
         I = I-1
         J = 2*I-1
    end
    if J > 1
        retval = 0.0
    end
    return retval
end

function PSI(M::Integer, K::Integer)
    FACTR1=0.0
    FACTR2=0.0
    N=2*K+M
    J=1
    if K>=1
        while J <= K
            FACTR1 = FACTR1+1.0/J
            J += 1
        end
    end
    I = 1
    while I <= N
        FACTR2 = FACTR2 + 1.0/(2.0*I-1.0)
        I=I+1
    end
    if N != 0
        FACTR2 = 0.0
    end

    return FACTR1-2.0*FACTR2

end



function FRAC(Z::Real, M::Integer, NMAX::Integer, ε::Real, TINYSQ::Real)
    N = NMAX
    MM = 0
    DZ2 = Z + Z
    DN0 = N + M
    DN1 = N + 1.0
    DN2 = DN0 + 0.5
    DN3 = DN0 - 0.5
    DN4 = 1.0 - M - M
    B = 2.0*DN1*Z/DN2
    A = 1.0
    FC = TINYSQ
    C0 = FC
    D0 = 0.0
    while true
        D0 = B + A*D0
        if abs(D0)  < TINYSQ
            D0 = TINYSQ
        end
        C0 = B + A/C0
        if abs(C0) < TINYSQ
            C0 = TINYSQ
        end
        D0 = 1.0/D0
        DELTA = C0*D0
        FC = FC*DELTA
        MM = MM + 1
        A = -(1.0 + DN4/(DN3 + MM))
        B = DZ2*(DN1 + MM)/(DN2 + MM)
        if MM < 1000000
            if abs(DELTA-1.0) > ε
                continue
            else
                return FC
            end
        elseif MM == 1000000
            println("CF CONVERGENCE FAILS")
            return FC
        end
    end
    return FC
end

function FACTOR(M::Integer, K::Integer)
    retval = 1.0
    if K >= 1
        X1 = M + 0.5
        N = 2*K - 1
        I = K
        J = N
        while J >= 0
            retval = retval*(J+X1)/I/I
            J = J-1
            I = I -1
            if I == 0
                I = 1
            end
        end
    end
    return retval
end


function EXPAN(Z::Real, MODE::Integer, IPRE::Integer, OVER::Real, QARGU::Real, M::Integer)
    PRECI = [1.e-13 1.e-9]
    πSQ=sqrt(π)
    DB = 2.0*log(2.0)
    FL = M/2.0
    CC = abs(trunc(FL)-FL)
    if (CC < 0.4)
        AR = 1.0
    else
        AR = -1.0
    end
    DZ = 1.0
    for I = 1:M
        DZ = DZ/QARGU
    end
    if MODE == 0
        GAMMA = GAMMAH(M, OVER)*AR*πSQ
    else
        GAMMA = AR
    end
    DFAC = 2.0/π*DZ*GAMMA/πSQ
    DF1 = log(2.0*Z)
    A0 = 1.0/sqrt(2.0*Z)
    Z2I = 1.0/(Z*Z)
    DELTA = 1.0
    SUM = 0.0
    K = 0
    DA2 = FACTOR(M, 0)
    DA1 = DB + PSI(M, 0)
    while true
        DELTA = (DF1+DA1)*DA2*A0
        SUM = SUM+DELTA
        DCC = 0.5+M+2.0*K
        DCCP = DCC + 1.0
        DKK = K + 1.0
        DA2 = DA2*DCCP*DCC/(DKK*DKK)
        DA1 = DA1+1.0/DKK-1.0/DCCP-1.0/DCC
        K = K + 1
        A0 = A0*0.25*Z2I
        if abs(DELTA/SUM) > PRECI[IPRE]
            continue
        else
            break
        end
    end
    return SUM*DFAC
end


function FRACPS(QZ::Real, M::Integer, N::Integer, ε::Real, TINYSQ::Real)
    MM = 0
    DN2 = N*N
    DZ2 = QZ+QZ
    DM = M-0.5
    B = DZ2*M
    A = DN2-DM*DM
    FC = TINYSQ
    C0 = FC
    D0 = 0.0
    while true
        D0 = B + A*D0
        if abs(D0) < TINYSQ
            D0 = TINYSQ
        end
        C0 = B + A/C0
        if abs(C0) < TINYSQ
            C0 = TINYSQ
        end
        D0 = 1.0/D0
        DELTA = C0*D0
        FC = FC*DELTA
        MM = MM + 1
        A = DN2-(MM+DM)*(MM+DM)
        B = DZ2*(MM+M)
        if MM < 1000000
            if  abs(DELTA-1.0) > ε
                continue
            else
                return FC
            end
        elseif MM == 1000000
            println("CF CONVERGENCE FAILS")
            return FC
        end
    end
    return FC
end

function FACTCO(N::Integer, PL::Real, M::Integer)
    retval = 1.0/PL
    X1 =  M + 0.50
    X2 = -M + 0.50
    J = N
    while J >= 0
        retval = retval*(J+X1)/(J+X2)
        J -= 1
    end
    return retval
end

function ELLIP2(AK)
    Q=(1.0-AK)*(1.0+AK)
    return DRF(Q)-(AK)^2*DRD(Q)/3.0
end



function ELLIP1(AK)
    return DRF((1.0-AK)*(1.0+AK))
end

function DRF(Y::Real)
    ERRTOL = 1.e-8
    C1 = 1.0/24.0
    C2 = 3.0/44.0
    C3 = 1.0/14.0

    DRF = 0.0
    XN = 0.
    YN = Y
    ZN = 1.0

    MU = 0.0
    XNDEV = 0.0
    YNDEV = 0.0
    ZNDEV = 0.0

    while true
        MU = (XN + YN + ZN)/3.0
        XNDEV = 2.0 - (MU+XN)/MU
        YNDEV = 2.0 - (MU+YN)/MU
        ZNDEV = 2.0 - (MU+ZN)/MU
        EPSLON = max(abs(XNDEV), abs(YNDEV), abs(ZNDEV))
        if EPSLON < ERRTOL
            break
        end
        XNROOT = sqrt(XN)
        YNROOT = sqrt(YN)
        ZNROOT = sqrt(ZN)
        LAMDA = XNROOT*(YNROOT + ZNROOT) + YNROOT*ZNROOT
        XN = (XN + LAMDA)*0.250
        YN = (YN + LAMDA)*0.250
        ZN = (ZN + LAMDA)*0.250
    end
    E2 = XNDEV*YNDEV - ZNDEV*ZNDEV
    E3 = XNDEV*YNDEV*ZNDEV
    S  = 1.0 + (C1*E2-0.10-C2*E3)*E2 + C3*E3
    return S/sqrt(MU)
end

function DRD(Y::Real)
    ERRTOL = 1.e-8
    C1 = 3.0/14.0
    C2 = 1.0/6.0
    C3 = 9.0/22.0
    C4 = 3.0/26.0

    DRD = 0.0

    XN = 0.0
    YN = Y
    ZN = 1.0
    SIGMA = 0.0
    POWER4 = 1.0

    MU = 0.0
    XNDEV = 0.0
    YNDEV = 0.0
    ZNDEV = 0.0

    while true
        MU = (XN + YN + 3.0*ZN)*0.20
        XNDEV = (MU-XN)/MU
        YNDEV = (MU-YN)/MU
        ZNDEV = (MU-ZN)/MU
        EPSLON = max(abs(XNDEV), abs(YNDEV), abs(ZNDEV))
        if EPSLON < ERRTOL
            break
        end
        XNROOT = sqrt(XN)
        YNROOT = sqrt(YN)
        ZNROOT = sqrt(ZN)
        LAMDA = XNROOT*(YNROOT + ZNROOT) + YNROOT*ZNROOT
        SIGMA = SIGMA + POWER4/(ZNROOT*(ZN+LAMDA))
        POWER4 = POWER4*0.250
        XN = (XN+LAMDA)*0.250
        YN = (YN+LAMDA)*0.250
        ZN = (ZN+LAMDA)*0.250
    end

    EA = XNDEV*YNDEV
    EB = ZNDEV*ZNDEV
    EC = EA - EB
    ED = EA - 6.0*EB
    EF = ED + EC + EC
    S1 = ED*(-C1+0.250*C3*ED-1.50*C4*ZNDEV*EF)
    S2 = ZNDEV*(C2*EF+ZNDEV*(-C3*EC+ZNDEV*C4*EA))
    return 3.0*SIGMA + POWER4*(1.0 + S1 + S2)/(MU*sqrt(MU))
 end 
