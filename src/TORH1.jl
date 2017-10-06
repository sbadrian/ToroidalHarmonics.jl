
"""
    TORH1(z, m, nmax; ε = 1e-14, ipre = 1, mode = 0)

Computes the toroidal harmonics P_{n-1/2}^m(z) and Q_{n-1/2}^m(z) of order m
from the lowest degree n = 0 to n = nmax
(Legendre polynomials with degree  n = -1/2 to n = -1/2 + nmax).

# Arguments
- ε = 1e-14: the desired precision
- ipre = 1: required precision in the evaluation of the toroidal harmonics.
        if ipre = 1, precision=10**(-12) (taking ε<10**(-12))
        if ipre = 2, precision=10**(-8)  (taking ε<10**(-8))
 - isgamma  = true: enforces a multiplication with with gamma(M + 1/2)
 - mode = 0: can also assume 1 and 2. (For mode = 1,
            results are divided by gamma(m + 1/2), for mode = 2 is good for
            high m's and for z > 20. )
"""
function TORH1(Z::Real, M::Integer, nmax::Integer;
                ε = 1e-14, ipre = 1, mode = 0, isgamma=true)

    const tiny = realmin(typeof(Z))*10.0^20 #1e-290

    ####################################################################
    # The dimension of QLMM (internal array) must be greater than m
    ####################################################################
    QLMM  = zeros(10000 + M*2)
    PL = zeros(nmax+2)
    QL =  zeros(nmax+2)

    over = 1/tiny
    tinySQ = sqrt(tiny)
    if (ipre != 1) && (ipre != 2)
        error("ipre must be 1 or 2")
    end

    PR = [0.22 0.12]
    nmaxP = nmax
    ####################################################################
    #   ε: required accuracy for the continued fraction
    #      (modified Lentz)
    #   tiny: small parameter to prevent overflows in the CF
    #         (close to the underflow limit)
    ####################################################################
    if Z <= 1.0
        error("Improper argument. Z must be greater than 1")
    end

    QZ=Z
    πSQ=sqrt(π)
    DPπ=sqrt(2.0)/πSQ

    FL=M/2.0
    CC=abs(trunc(FL)-FL)
    if CC < 0.4
      AR = 1.0
    else
      AR = -1.0
    end

    ####################################################################
    ##  wh choose expansion or cf for PL[0] depending on the values
    ##  of Z,M and mode
    ####################################################################
    ICAL = 1
    if M != 0
        if Z/M > PR[ipre]
            ICAL = 2
        end

        if Z < 5.0
            ICAL = 1
        elseif Z > 20.0
            if (mode != 2) && (ICAL == 1)
                ICAL = 0
            end
        end
        if ICAL == 0
            error("You must choose mode = 2")
        end
    else
        if Z  < 5.0
            ICAL = 1
        else
            ICAL = 2
        end
    end
    ####################################################################                                                                       C
    #   we use the code if nmax is greater than or equal to 2
    ####################################################################
    # if nmaxP < 2
    #    nmaxP = 2
    # end

    if mode == 0
        GAMMA = GAMMAH(M, over)*AR*πSQ
        if abs(GAMMA) < tiny
            error("M is too large for mode = 0, better try mode = 1")
        end
    else
        GAMMA = AR
    end
    ####################################################################
    #      We evalue the continued fraction using
    #      Lentz-Thompson
    ####################################################################
    FC = FRAC(Z, M, 0, ε, tinySQ)
    QDC1 = QZ*QZ-1.0
    QARGU = QZ/sqrt(QDC1)
    DFAC1 = DPπ*GAMMA/π
    DFAC2 = GAMMA/DPπ

    ####################################################################
    ##   We  evaluate Q_{-1/2}, Q^{1}_{-1/2}
    ##   using SLATEC routines for elliptic functions
    ####################################################################
    ARGU1 = sqrt(2.0/(Z+1.0))
    QLMM[1]=ARGU1*ELLIP1(ARGU1)
    QLMM[2]=-1.0/sqrt(2.0*(QZ-1.0))*ELLIP2(ARGU1)

    ####################################################################
    ## we apply forward recurrence in m for Q'S
    ####################################################################
    MP = 1
    if mode == 0
        while (MP <= M) && (abs(QLMM[MP+1]) < over)
            QLMM[MP + 2] = -2.0*MP*QARGU*QLMM[MP+1]-(MP-0.5)*(MP-0.5)*QLMM[MP]
            MP += 1
        end
        if (MP - 1) != M
           error("M is too large for mode = 0, better try mode = 1")
        end
    else
        QLMM[1] = QLMM[1]/πSQ
        QLMM[2] = QLMM[2]*2.0/πSQ
        while (MP <= M) && (abs(QLMM[MP + 1]) < over)
            D1 = MP + 0.5
            QLMM[MP + 2] = -2.0*MP*QARGU*QLMM[MP + 1]/D1 - (MP-0.5)*QLMM[MP]/D1
            MP = MP + 1
        end
        if (MP - 1) != M
            error("M is too large for mode = 1, 2")
        end
    end
    nmmax = M
    DFACQS = -(GAMMA/QLMM[M + 2])*GAMMA/π

    ####################################################################
    ##  Evaluation of PL[0]
    ####################################################################

    if ICAL == 1
        ####################################################################
        ##   We calculate the CF for P's
        ####################################################################
        FCP = FRACPS(QARGU, M + 1, 0, ε, tinySQ)
        DD = M + 0.5
        if mode != 0
            FCP = FCP/DD
            DFACQS = DFACQS/DD
        end
        PL[1] = DFACQS/sqrt(QDC1)/(1.0-FCP*QLMM[M + 1]/QLMM[M + 2])
    else
        PL0 = EXPAN(Z, mode, ipre, over, QARGU, M)
        PL[1] = PL0
    end
    QM0 = QLMM[M + 1]
    DFAC3 = (GAMMA/QM0)*GAMMA/π/(0.5 - M)
    PL[2] = PL[1]*FC + DFAC3
    NP = 1
    ####################################################################                                                                         C
    #   We use the current relations for the P's
    ####################################################################
    while (NP <= nmaxP) && (abs((NP - M + 0.5)*PL[NP + 1]) < over)
        PL[NP + 2] = (2.0*NP*Z*PL[NP + 1]-(NP + M - 0.5)*PL[NP])/(NP - M + 0.5)
        NP += 1
    end

    nmaxP = NP - 1
    newn = nmaxP
    DFAC4 = (FACTCO(nmaxP,PL[nmaxP + 2], M)*GAMMA)*GAMMA/π/(nmaxP + M + 0.5)

    ####################################################################
    #  We evalulate the continued fraction using Lentz-Thompson                                                                                  C
    ####################################################################
    FC = FRAC(Z, M, nmaxP, ε, tinySQ)

    ####################################################################
    #     Evaluation of PL[nmax + 2] and PL[nmax + 1] using
    #     the Wronskian W{PL[nmax + 1], QL[nmax + 1]},
    #     the known values of PL[nmax + 2] and PL[nmax + 1]
    #     the value of H = QL[nmax + 2]/QL[nmax + 1]
    ####################################################################
    QL[nmaxP + 1] = DFAC4/(1.0-FC*PL[nmaxP + 1]/PL[nmaxP + 2])
    QL[nmaxP + 2] = QL[nmaxP + 1]*FC

    ####################################################################
    #   We use the backward recurrence relation
    ####################################################################
    for I = 1:nmaxP
        NP = nmaxP + 1 - I
        N = NP - 1
        QL[N+1]=((NP+NP)*Z*QL[NP + 1]-(NP-M+0.5)*QL[NP + 2])/(NP + M - 0.5)
    end

    if nmax == 0
        return PL[1:1], QL[1:1]
    else
        return PL, QL
    end

end
