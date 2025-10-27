PROGRAM dvl5                       IMPLICIT NONE
  !----------------------------- grid/config -----------------------------
  INTEGER, PARAMETER :: NZ = 60            ! vertical levels (surface = 1)
  DOUBLE PRECISION, PARAMETER :: DZ = 0.5D0 ! layer thickness [m]   DOUBLE PRECISION, PARAMETER :: HMAX = NZ*DZ
                                   !----------------------------- time setup ------------------------------
  INTEGER, PARAMETER :: NDAYS = 3  DOUBLE PRECISION, PARAMETER :: DT = 300.D0         ! timestep [s] (5 min)
  INTEGER, PARAMETER :: STEPS_PER_DAY = INT(86400.D0/DT + 0.5D0)
  INTEGER, PARAMETER :: NSTEPS = NDAYS*STEPS_PER_DAY
                                   !--------------------------- physical params ---------------------------
  DOUBLE PRECISION, PARAMETER :: FK = -2.D0*7.2921159D-5*DSIN(-37.D0*3.14159265358979D0/180.D0)
  !   fk ~ f = 2 Ω sin φ, φ≈ -37° (Melbourne); sign gives clockwise inertial rot in SH               DOUBLE PRECISION, PARAMETER :: G = 9.81D0
  DOUBLE PRECISION, PARAMETER :: RHO0 = 1025.D0
  DOUBLE PRECISION, PARAMETER :: ALPHA_T = 2.0D-4 ! thermal expansivity [1/K]                        DOUBLE PRECISION, PARAMETER :: CPW = 3985.D0    ! seawater cp [J/(kg K)]
                                   ! Shortwave penetration: two exponential scales
  DOUBLE PRECISION, PARAMETER :: F1 = 0.6D0       ! fraction in first e-fold                         DOUBLE PRECISION, PARAMETER :: ZETA1 = 0.35D0   ! m, very near-surface
  DOUBLE PRECISION, PARAMETER :: ZETA2 = 23.0D0   ! m, blue/green band

  ! Background & mixing
  DOUBLE PRECISION, PARAMETER :: KM0 = 5.D-5      ! background K [m^2/s]
  DOUBLE PRECISION, PARAMETER :: KH0 = 5.D-5                        DOUBLE PRECISION, PARAMETER :: RIC = 0.25D0     ! critical Richardson
                                   !----------------------------- arrays/vars ------------------------------
  DOUBLE PRECISION, ALLOCATABLE :: ZH(:), ZI(:)    ! layer centers & interfaces
  DOUBLE PRECISION, ALLOCATABLE :: U(:), V(:), T(:), N2(:)          DOUBLE PRECISION, ALLOCATABLE :: KM(:), KH(:), QT(:) ! mixing & heat source
  DOUBLE PRECISION :: TIME, HOD, QSW, TAU, ZML, KE_TEND, WW, SP
  DOUBLE PRECISION :: UN, VN, DU, DV, S2, RICH, FAC
  INTEGER :: K, n, DAYIDX
  CHARACTER(LEN=64) :: FNAME
                                   !----------------------------- init grids ------------------------------                           ALLOCATE(ZH(NZ), ZI(NZ+1))
  DO K = 1, NZ
     ZH(K) = (K-0.5D0)*DZ              ! depth positive downward [m]
  END DO                           DO K = 1, NZ+1                      ZI(K) = (K-1)*DZ
  END DO

  !----------------------------- init state ------------------------------
  ALLOCATE(U(NZ), V(NZ), T(NZ), N2(NZ), KM(NZ), KH(NZ), QT(NZ))     U = 0.D0; V = 0.D0
  T = 20.0D0 - 0.02D0*ZH             ! weakly warmer at top (stable)
  CALL compute_n2(T, ZH, N2, ALPHA_T, G, DZ)                        KM = KM0; KH = KH0               QT = 0.D0

  ! Diagnostics init               ZML = 10.D0
  KE_TEND = 0.D0; WW = 0.D0; SP = 0.D0
                                   !----------------------------- time loop -------------------------------                           TIME = 0.D0
  DO n = 1, NSTEPS
     HOD     = MOD(TIME/3600.D0, 24.D0)        ! hour of day [0..24)
     DAYIDX  = INT(TIME/86400.D0) + 1

     ! Surface forcing (synthetic): clear-sky SW, weak wind stress     QSW = diurnal_qsw(HOD)                     ! W/m^2, >=0 during day                                 TAU = 0.04D0                               ! N/m^2 (weak wind)
                                      ! Build penetrative heating profile [K/s] in QT
     CALL shortwave_heating(QSW, F1, ZETA1, ZETA2, ZH, DZ, RHO0, CPW, QT)

     ! Shear & mixing coefficients via gradient Richardson limiter
     DO K = 2, NZ
        DU   = U(K) - U(K-1)             DV   = V(K) - V(K-1)
        S2   = (DU*DU + DV*DV)/(DZ*DZ) + 1.D-10                           RICH = MAX(0.D0, 0.5D0*(N2(K)+N2(K-1)) / S2)
        FAC  = 1.D0/(1.D0 + 5.D0*(RICH/RIC))   ! simple limiter
        KM(K) = KM0 + 0.03D0*FAC               ! toy dependence           KH(K) = KH0 + 0.05D0*FAC
     END DO
     KM(1) = KM(2); KH(1) = KH(2)
     ! Apply surface shear (wind stress) as a boundary momentum flux                                    U(1) = U(1) + (TAU/RHO0)*DT/DZ
     ! (Assume stress aligned with +u; set V stress to 0 here)
     ! Optionally: V(1) = V(1) + (TAU_V/RHO0)*DT/DZ
                                      ! Coriolis predictor for (U,V)
     DO K = 1, NZ                        UN = U(K) - FK * V(K) * DT
        VN = V(K) + FK * U(K) * DT                                        U(K) = UN
        V(K) = VN                     END DO
                                      ! Vertical diffusion step (explicit FTCS with harmonic KM/KH)     CALL diffuse_tridiag(U, KM, DZ, DT)                               CALL diffuse_tridiag(V, KM, DZ, DT)
     CALL diffuse_tridiag_scalar(T, KH, DZ, DT)                                                         ! Add penetrative heating
     T = T + QT*DT                                                     ! Update buoyancy frequency
     CALL compute_n2(T, ZH, N2, ALPHA_T, G, DZ)                   
     ! Mixed layer depth (simple: ΔT from surface > 0.2 K)             ZML = mld_deltaT(T, ZH, 0.2D0)

     ! Crude diagnostics: kinetic energy tendency proxy, shear prod, …                                  KE_TEND = 0.D0
     DO K = 1, NZ                        KE_TEND = KE_TEND + 0.5D0*(U(K)*U(K)+V(K)*V(K))*DZ             END DO
     SP = 0.D0                        DO K = 2, NZ                        DU = U(K)-U(K-1); DV = V(K)-V(K-1)                                SP = SP + KM(K)*( (DU/DZ)**2 + (DV/DZ)**2 )*DZ
     END DO
     WW = 0.D0   ! placeholder (no conv. vert vel in this toy)
                                      ! Print hourly
     IF (MOD(n, INT(3600.D0/DT)) == 0) THEN
        WRITE(*,'(A,I0,2X,A,F6.2,2X,A,F6.2,2X,A,F7.3,2X,A,F6.2,2X,A,F7.4,2X,A,F7.4,2X,A,F8.3)') &            'Day', DAYIDX, 'Hr', HOD, 'MLD[m]', ZML, 'U_sfc[m/s]', &
          SQRT(U(1)**2 + V(1)**2), 'dT0[K]', MAX(0.D0, T(1)-T(2)), 'KEtend', KE_TEND, 'WW', WW, 'SP', SP
     END IF                                                            ! Dump simple profile file at 18:00 each day
     IF (ABS(HOD-18.D0) .LT. 1.D-6) THEN                                  WRITE(FNAME,'(A,I0,A,I0,A)') 'dwl_TU_day', DAYIDX, '_', INT(HOD), 'h.dat'
        CALL write_profile(FNAME, ZH, U, V, T)                         END IF                      
     TIME = TIME + DT
  END DO                         
  ! Final dump                     CALL write_profile('dwl_final.dat', ZH, U, V, T)

CONTAINS                         
  SUBROUTINE compute_n2(T, ZH, N2, ALPHA_T, G, DZ)
    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(IN)  :: T(:), ZH(:), ALPHA_T, G, DZ      DOUBLE PRECISION, INTENT(OUT) :: N2(SIZE(T))                      INTEGER :: k, nz
    DOUBLE PRECISION :: dTdz         nz = SIZE(T)                     N2(1) = 0.D0                     DO k = 2, nz
       dTdz = (T(k)-T(k-1))/DZ          N2(k) = -G * (-ALPHA_T) * dTdz   ! ρ ≈ ρ0 (1 - α_T ΔT); N2 = -(g/ρ) dρ/dz
    END DO                         END SUBROUTINE compute_n2      
  SUBROUTINE shortwave_heating(QSW, F1, ZETA1, ZETA2, ZH, DZ, RHO0, CPW, QT)
    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(IN)  :: QSW, F1, ZETA1, ZETA2, ZH(:), DZ, RHO0, CPW                       DOUBLE PRECISION, INTENT(OUT) :: QT(:)   ! [K/s]                  INTEGER :: k, nz                 DOUBLE PRECISION :: Fz1a, Fz1b, Fz2a, Fz2b, dF
    nz = SIZE(ZH)
    DO k = 1, nz                        ! Layer-integrated absorbed SW over [z1,z2] divided by (ρ cp Δz)                                   Fz1a = F1 * DEXP( -(ZH(k)-0.5D0*DZ)/ZETA1 )
       Fz1b = F1 * DEXP( -(ZH(k)+0.5D0*DZ)/ZETA1 )
       Fz2a = (1.D0-F1) * DEXP( -(ZH(k)-0.5D0*DZ)/ZETA2 )                Fz2b = (1.D0-F1) * DEXP( -(ZH(k)+0.5D0*DZ)/ZETA2 )                dF   = (Fz1a - Fz1b) + (Fz2a - Fz2b)   ! fraction absorbed in layer
       QT(k) = (QSW * dF) / (RHO0 * CPW * 1.D0)   ! W/m^2 * frac / (ρ cp)  -> K/s (since per m^2 over 1 m^2 column)
       QT(k) = QT(k) / DZ                         ! distribute per thickness -> K/s per layer-avg      END DO
  END SUBROUTINE shortwave_heating

  FUNCTION diurnal_qsw(hod) RESULT(q)                                 IMPLICIT NONE
    DOUBLE PRECISION, INTENT(IN) :: hod                               DOUBLE PRECISION :: q, mu
    ! Clear-sky toy: sunrise 6, sunset 18, peak ~ 800 W/m^2 at noon                                    IF (hod >= 6.D0 .AND. hod <= 18.D0) THEN
       mu = DSIN( (hod-6.D0)*3.14159265358979D0/12.D0 )   ! 0..1..0
       q  = 800.D0 * mu              ELSE                                q = 0.D0
    END IF
  END FUNCTION diurnal_qsw       
  SUBROUTINE diffuse_tridiag(X, Kcoef, DZ, DT)                        IMPLICIT NONE
    DOUBLE PRECISION, INTENT(INOUT) :: X(:)                           DOUBLE PRECISION, INTENT(IN)    :: Kcoef(:), DZ, DT
    INTEGER :: nz, k                 DOUBLE PRECISION, ALLOCATABLE :: a(:), b(:), c(:), rhs(:), xm(:)
    DOUBLE PRECISION :: kmh, kph, r
    nz = SIZE(X)                     ALLOCATE(a(nz), b(nz), c(nz), rhs(nz), xm(nz))                    rhs = X

    ! Build tri-diagonal for implicit diffusion: dX/dt = d/dz (K dX/dz)                                ! Crank–Nicolson (theta=0.5)
    a = 0.D0; b = 1.D0; c = 0.D0     DO k = 2, nz-1
       kmh = 0.5D0*(Kcoef(k-1)+Kcoef(k))                                 kph = 0.5D0*(Kcoef(k)+Kcoef(k+1))
       r = 0.5D0*DT/(DZ*DZ)
       a(k) = -r*kmh                    b(k) = 1.D0 + r*(kmh+kph)
       c(k) = -r*kph
    END DO                       
    ! Boundary conditions: no-flux at bottom, flux specified at surface already applied in X           b(1)   = 1.D0 + (0.5D0*DT/(DZ*DZ))*(0.5D0*(Kcoef(1)+Kcoef(2)))    c(1)   = - (0.5D0*DT/(DZ*DZ))*(0.5D0*(Kcoef(1)+Kcoef(2)))
    a(nz)  = - (0.5D0*DT/(DZ*DZ))*(0.5D0*(Kcoef(nz-1)+Kcoef(nz)))     b(nz)  = 1.D0 + (0.5D0*DT/(DZ*DZ))*(0.5D0*(Kcoef(nz-1)+Kcoef(nz)))                             
    CALL solve_tridiag(a, b, c, rhs, xm)
    X = xm                           DEALLOCATE(a,b,c,rhs,xm)       END SUBROUTINE diffuse_tridiag
                                   SUBROUTINE diffuse_tridiag_scalar(T, Kcoef, DZ, DT)
    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(INOUT) :: T(:)                           DOUBLE PRECISION, INTENT(IN)    :: Kcoef(:), DZ, DT               CALL diffuse_tridiag(T, Kcoef, DZ, DT)
  END SUBROUTINE diffuse_tridiag_scalar                           
  SUBROUTINE solve_tridiag(a,b,c,d,x)                                 IMPLICIT NONE
    DOUBLE PRECISION, INTENT(IN)  :: a(:), b(:), c(:), d(:)
    DOUBLE PRECISION, INTENT(OUT) :: x(:)
    INTEGER :: n, i
    DOUBLE PRECISION, ALLOCATABLE :: cp(:), dp(:)
    n = SIZE(d)
    ALLOCATE(cp(n), dp(n))
    cp(1) = c(1)/b(1)                dp(1) = d(1)/b(1)
    DO i = 2, n                         cp(i) = c(i)/(b(i)-a(i)*cp(i-1))
       dp(i) = (d(i)-a(i)*dp(i-1))/(b(i)-a(i)*cp(i-1))
    END DO
    x(n) = dp(n)                     DO i = n-1, 1, -1
       x(i) = dp(i) - cp(i)*x(i+1)
    END DO                           DEALLOCATE(cp,dp)
  END SUBROUTINE solve_tridiag                                      FUNCTION mld_deltaT(T, ZH, thresh) RESULT(zml)                      IMPLICIT NONE
    DOUBLE PRECISION, INTENT(IN) :: T(:), ZH(:), thresh
    DOUBLE PRECISION :: zml          INTEGER :: k, nz
    DOUBLE PRECISION :: T0
    nz = SIZE(T)                     T0 = T(1)
    zml = ZH(nz)
    DO k = 1, nz                        IF (ABS(T(k)-T0) > thresh) THEN
          zml = ZH(k)
          EXIT                          END IF
    END DO                         END FUNCTION mld_deltaT

  SUBROUTINE write_profile(fname, ZH, U, V, T)
    IMPLICIT NONE                    CHARACTER(*), INTENT(IN) :: fname
    DOUBLE PRECISION, INTENT(IN) :: ZH(:), U(:), V(:), T(:)
    INTEGER :: k, nz, ios
    OPEN(10, FILE=fname, STATUS='REPLACE', ACTION='WRITE', IOSTAT=ios)
    IF (ios /= 0) RETURN
    WRITE(10,'(A)') '# z[m]  U[m/s]  V[m/s]  T[C]'
    nz = SIZE(ZH)                    DO k = 1, nz                        WRITE(10,'(F8.3,1X,F10.6,1X,F10.6,1X,F8.4)') ZH(k), U(k), V(k), T(k)
    END DO
    CLOSE(10)                      END SUBROUTINE write_profile   
END PROGRAM dvl5
