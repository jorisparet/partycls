MODULE compute

  IMPLICIT NONE

  ! Constant Pi
  REAL(8),  PARAMETER :: pi  = 4.0_8 * ATAN(1.0_8)

  ! Tabulated factorials
  INTEGER(16), PARAMETER, DIMENSION(33) :: factorial = (/ &
    1_16, &
    1_16, &
    2_16, &
    6_16, &
    24_16, &
    120_16, &
    720_16, &
    5040_16, &
    40320_16, &
    362880_16, &
    3628800_16, &
    39916800_16, &
    479001600_16, &
    6227020800_16, &
    87178291200_16, &
    1307674368000_16, &
    20922789888000_16, &
    355687428096000_16, &
    6402373705728000_16, &
    121645100408832000_16, &
    2432902008176640000_16, &
    51090942171709440000_16, &
    1124000727777607680000_16, &
    25852016738884976640000_16, &
    620448401733239439360000_16, &
    15511210043330985984000000_16, &
    403291461126605635584000000_16, &
    10888869450418352160768000000_16, &
    304888344611713860501504000000_16, &
    8841761993739701954543616000000_16, &
    265252859812191058636308480000000_16, &
    8222838654177922817725562880000000_16, &
    263130836933693530167218012160000000_16 &
  /)
  
CONTAINS

  !!!!!!!!!!! PBC !!!!!!!!!!!
  SUBROUTINE pbc(d_ij, box, hbox)
      REAL(8), INTENT(inout)  :: d_ij(:)
      REAL(8), INTENT(in)     :: box(:)
      REAL(8), INTENT(in)     :: hbox(:)
      WHERE (ABS(d_ij) > hbox)
          d_ij = d_ij - SIGN(box ,d_ij)
      END WHERE
  END SUBROUTINE pbc
  

  SUBROUTINE radial_histogram(pos_0, pos_1, idx_0, idx_1, box, r, dr, hist)
    ! Parameters
    REAL(8), INTENT(in)     :: pos_0(:,:), pos_1(:,:), box(:), r(:), dr
    INTEGER(8), INTENT(in)  :: idx_0(:), idx_1(:)
    INTEGER(8), INTENT(out) :: hist(SIZE(idx_0), SIZE(r))
    ! Variables
    INTEGER(8) :: i, j, idx0, idx1, np0, np1, bin, nbins
    REAL(8)    :: r_ij(SIZE(box)), d_ij_sq, hbox(SIZE(box)), rmin
    ! Computation
    np0 = SIZE(idx_0)
    np1 = SIZE(idx_1)
    hbox = box / 2.0
    nbins = SIZE(r)
    rmin = r(1)
    ! initialize to zero
    hist = 0
    ! iterate over main group
    DO i=1,np0
      idx0 = idx_0(i)
      ! iterate over second group
      DO j=1,np1
        idx1 = idx_1(j)
        ! do not consider particle `i` with itself
        IF (idx0 /= idx1) THEN
          r_ij(:) = pos_0(:,i) - pos_1(:,j)
          CALL pbc(r_ij, box, hbox)
          d_ij_sq = SUM(r_ij**2)
          ! histogram
          bin = FLOOR(SQRT(d_ij_sq) / dr) - FLOOR(rmin / dr) + 1
          IF (bin >= 1 .AND. bin <= nbins) hist(i, bin) = hist(i, bin) + 1
        END IF
      END DO
    END DO
  END SUBROUTINE

  !!!!!!!!!! ANGULAR HISTOGRAM (INDIVIDUAL) !!!!!!!!!!
  SUBROUTINE angular_histogram(idx_i, pos_i, pos_all, neigh_i, box, nbins, dtheta, hist)
    ! Parameters
    INTEGER(8), INTENT(in)  :: idx_i, neigh_i(:), nbins
    REAL(8), INTENT(in)     :: pos_i(:), pos_all(:,:), box(:), dtheta
    INTEGER(8), INTENT(out) :: hist(nbins)
    ! Variables
    INTEGER(8) :: j, k, idx_j, idx_k, n_neigh_i, bin
    REAL(8)    :: r_ij(SIZE(box)), r_ik(SIZE(box)), d_ij, d_ik, hbox(SIZE(box))
    REAL(8)    :: dotprod, prod, costheta, theta
    ! Computation
    hbox = box / 2.0
    n_neigh_i = SIZE(neigh_i)
    ! initialize to zero
    hist = 0
    ! First neighbor: j
    DO j=1,n_neigh_i
      idx_j = neigh_i(j) + 1 ! python index shift
      IF (idx_j /= idx_i+1) THEN ! pass if j=i
        r_ij(:) = pos_i - pos_all(:,idx_j)
        CALL pbc(r_ij, box, hbox)
        d_ij = SQRT(SUM(r_ij**2))
        ! Second neighbor: k
        DO k=1,n_neigh_i
          idx_k = neigh_i(k) + 1 ! python index shift
          IF (idx_k /= idx_i+1 .AND. idx_k /= idx_j) THEN ! pass if k=i or k=j
            r_ik(:) = pos_i - pos_all(:,idx_k)
            CALL pbc(r_ik, box, hbox)
            d_ik = SQRT(SUM(r_ik**2))
            ! angle (k,i,j)
            dotprod = SUM(r_ij*r_ik)
            prod = d_ij*d_ik
            costheta = dotprod/prod
            ! enforce cos(theta) >= -1
            IF (costheta <= 0.0) THEN
              costheta = DMAX1(-1.0_8,costheta)
            END IF
            ! enforce cos(theta) <= 1
            IF (costheta > 0.0) THEN
              costheta = DMIN1(1.0_8,costheta)
            END IF
            theta = ACOS(costheta)*180.0/pi
            ! binning
            bin = FLOOR( theta/dtheta ) + 1
            IF (bin <= nbins) hist(bin) = hist(bin) + 1
          END IF
        END DO
      END IF
    END DO
  END SUBROUTINE angular_histogram


  !!!!!!!!!! ANGULAR HISTOGRAM (ALL) !!!!!!!!!!
  SUBROUTINE angular_histogram_all(idx_0, pos_0, pos_all, neigh, neigh_number, &
                                   box, nbins, dtheta, hist)
    ! Parameters
    INTEGER(8), INTENT(in)  :: idx_0(:), neigh(:,:), neigh_number(:), nbins
    REAL(8), INTENT(in)     :: pos_0(:,:), pos_all(:,:), box(:), dtheta
    INTEGER(8), INTENT(out) :: hist(SIZE(idx_0), nbins)
    ! Variables
    INTEGER(8) :: i, j, k, idx_i, idx_j, idx_k, nn_i, bin
    REAL(8)    :: r_ij(SIZE(box)), r_ik(SIZE(box)), d_ij, d_ik, hbox(SIZE(box))
    REAL(8)    :: dotprod, prod, costheta, theta
    ! Computation
    hbox = box / 2.0
    hist = 0
    DO i=1,SIZE(idx_0)
      idx_i = idx_0(i) + 1 ! python index shift
      nn_i = neigh_number(i)
      ! first neighbor: j
      DO j=1,nn_i
        idx_j = neigh(i,j) + 1 ! python index shift
        IF (idx_j /= idx_i) THEN ! pass if j=i
          r_ij(:) = pos_0(:,i) - pos_all(:,idx_j)
          CALL pbc(r_ij, box, hbox)
          d_ij = SQRT(SUM(r_ij**2))
          ! second neighbor: k
          DO k=1,nn_i
            idx_k = neigh(i,k) + 1 ! python index shift
            IF (idx_k /= idx_i .AND. idx_k /= idx_j) THEN ! pass if k=i or k=j
              r_ik(:) = pos_0(:,i) - pos_all(:,idx_k)
              CALL pbc(r_ik, box, hbox)
              d_ik = SQRT(SUM(r_ik**2))
              ! angle (k,i,j)
              dotprod = SUM(r_ij*r_ik)
              prod = d_ij*d_ik
              costheta = dotprod/prod
              ! enforce cos(theta) >= -1
              IF (costheta <= 0.0) THEN
                costheta = DMAX1(-1.0_8,costheta)
              END IF
              ! enforce cos(theta) <= 1
              IF (costheta > 0.0) THEN
                costheta = DMIN1(1.0_8,costheta)
              END IF
              theta = ACOS(costheta)*180.0/pi
              ! binning
              bin = FLOOR( theta/dtheta ) + 1
              IF (bin <= nbins) THEN
                hist(i,bin) = hist(i,bin) + 1
              END IF
            END IF
          END DO
        END IF
      END DO
    END DO
  END SUBROUTINE angular_histogram_all

  
  SUBROUTINE smoothed_angular_histogram(idx_i, pos_i, pos_all, spe_i, spe_all, neigh_i, &
                                        cutoffs, pow, box, nbins, dtheta, hist)
    ! Parameters
    INTEGER(8), INTENT(in)  :: idx_i, spe_i, spe_all(:), neigh_i(:), pow, nbins
    REAL(8), INTENT(in)     :: pos_i(:), pos_all(:,:), cutoffs(:,:), box(:), dtheta
    REAL(8), INTENT(out)    :: hist(nbins)
    ! Variables
    INTEGER(8) :: j, k, idx_j, idx_k, n_neigh_i, bin
    REAL(8)    :: r_ij(SIZE(box)), r_ik(SIZE(box)), d_ij, d_ik, rc_ij, rc_ik, hbox(SIZE(box)), w_i
    REAL(8)    :: dotprod, prod, costheta, theta
    ! Computation
    hbox = box / 2.0
    n_neigh_i = SIZE(neigh_i)
    ! initialize to zero
    hist = 0.0
    ! First neighbor: j
    DO j=1,n_neigh_i
      idx_j = neigh_i(j) + 1 ! python index shift
      IF (idx_j /= idx_i+1) THEN ! pass if j=i
        r_ij(:) = pos_i - pos_all(:,idx_j)
        CALL pbc(r_ij, box, hbox)
        d_ij = SQRT(SUM(r_ij**2))
        ! Second neighbor: k
        DO k=1,n_neigh_i
          idx_k = neigh_i(k) + 1 ! python index shift
          IF (idx_k /= idx_i+1 .AND. idx_k /= idx_j) THEN ! pass if k=i or k=j
            r_ik(:) = pos_i - pos_all(:,idx_k)
            CALL pbc(r_ik, box, hbox)
            d_ik = SQRT(SUM(r_ik**2))
            ! angle (k,i,j)
            dotprod = SUM(r_ij*r_ik)
            prod = d_ij*d_ik
            costheta = dotprod/prod
            ! enforce cos(theta) >= -1
            IF (costheta <= 0.0) THEN
              costheta = DMAX1(-1.0_8,costheta)
            END IF
            ! enforce cos(theta) <= 1
            IF (costheta > 0.0) THEN
              costheta = DMIN1(1.0_8,costheta)
            END IF
            theta = ACOS(costheta)*180.0/pi
            ! binning
            bin = FLOOR( theta/dtheta ) + 1
            IF (bin <= nbins) THEN
              ! weights
              rc_ij = cutoffs(spe_i, spe_all(idx_j))
              rc_ik = cutoffs(spe_i, spe_all(idx_k))
              w_i = EXP( -( (d_ij/rc_ij)**pow + (d_ik/rc_ik)**pow ) )
              hist(bin) = hist(bin) + w_i
            END IF
          END IF
        END DO
      END IF
    END DO
  END SUBROUTINE smoothed_angular_histogram


  SUBROUTINE smoothed_angular_histogram_all(idx_0, pos_0, pos_all, spe_0, spe_all, & 
                                           neigh, neigh_number, cutoffs, pow, box, &
                                           nbins, dtheta, hist)
    ! Parameters
    INTEGER(8), INTENT(in)  :: idx_0(:), spe_0(:), spe_all(:), neigh(:,:), neigh_number(:)
    INTEGER(8), INTENT(in)  :: pow, nbins
    REAL(8), INTENT(in)     :: pos_0(:,:), pos_all(:,:), cutoffs(:,:), box(:), dtheta
    REAL(8), INTENT(out)    :: hist(SIZE(idx_0), nbins)
    ! Variables
    INTEGER(8) :: i, j, k, idx_i, idx_j, idx_k, nn_i, bin
    REAL(8)    :: r_ij(SIZE(box)), r_ik(SIZE(box)), d_ij, d_ik, rc_ij, rc_ik, hbox(SIZE(box)), w_i
    REAL(8)    :: dotprod, prod, costheta, theta
    ! Computation
    hbox = box / 2.0
    hist = 0.0
    DO i=1,SIZE(idx_0)
      idx_i = idx_0(i) + 1 ! python index shift
      nn_i = neigh_number(i)
      ! first neighbor: j
      DO j=1,nn_i
        idx_j = neigh(i,j) + 1 ! python index shift
        IF (idx_j /= idx_i) THEN ! pass if j=i
          r_ij(:) = pos_0(:,i) - pos_all(:,idx_j)
          CALL pbc(r_ij, box, hbox)
          d_ij = SQRT(SUM(r_ij**2))
          ! Second neighbor: k
          DO k=1,nn_i
            idx_k = neigh(i,k) + 1 ! python index shift
            IF (idx_k /= idx_i .AND. idx_k /= idx_j) THEN ! pass if k=i or k=j
              r_ik(:) = pos_0(:,i) - pos_all(:,idx_k)
              CALL pbc(r_ik, box, hbox)
              d_ik = SQRT(SUM(r_ik**2))
              ! angle (k,i,j)
              dotprod = SUM(r_ij*r_ik)
              prod = d_ij*d_ik
              costheta = dotprod/prod
              ! enforce cos(theta) >= -1
              IF (costheta <= 0.0) THEN
                costheta = DMAX1(-1.0_8,costheta)
              END IF
              ! enforce cos(theta) <= 1
              IF (costheta > 0.0) THEN
                costheta = DMIN1(1.0_8,costheta)
              END IF
              theta = ACOS(costheta)*180.0/pi
              ! binning
              bin = FLOOR( theta/dtheta ) + 1
              IF (bin <= nbins) THEN
                ! weights
                rc_ij = cutoffs(spe_0(i), spe_all(idx_j))
                rc_ik = cutoffs(spe_0(i), spe_all(idx_k))
                w_i = EXP( -( (d_ij/rc_ij)**pow + (d_ik/rc_ik)**pow ) )
                hist(i,bin) = hist(i,bin) + w_i
              END IF
            END IF
          END DO
        END IF
      END DO
    END DO
  END SUBROUTINE smoothed_angular_histogram_all
  
  !!!!!!!!!! TETRAHEDRALITY (INDIVIDUAL) !!!!!!!!!!
  SUBROUTINE tetrahedrality(idx_i, pos_i, pos_all, neigh_i, box, tetra)
    ! Parameters
    INTEGER(8), INTENT(in)  :: idx_i, neigh_i(:)
    REAL(8), INTENT(in)     :: pos_i(:), pos_all(:,:), box(:)
    REAL(8), INTENT(out)    :: tetra
    ! Variables
    INTEGER(8) :: j, k, idx_j, idx_k, n_neigh_i, N_ba
    REAL(8)    :: r_ij(SIZE(box)), r_ik(SIZE(box)), d_ij, d_ik, hbox(SIZE(box))
    REAL(8)    :: dotprod, prod, costheta_kij, costheta_tetra
    ! Computation
    hbox = box / 2.0
    n_neigh_i = SIZE(neigh_i)
    costheta_tetra = -0.333806859233771 ! cos(109.5°)
    N_ba = 0
    tetra = 0.0
    ! First neighbor: j
    DO j=1,n_neigh_i
      idx_j = neigh_i(j) + 1 ! python index shift
      IF (idx_j /= idx_i+1) THEN ! pass if j=i
        r_ij(:) = pos_i - pos_all(:,idx_j)
        CALL pbc(r_ij, box, hbox)
        d_ij = SQRT(SUM(r_ij**2))
        ! Second neighbor: k
        DO k=1,n_neigh_i
          idx_k = neigh_i(k) + 1 ! python index shift
          IF (idx_k /= idx_i+1 .AND. idx_k /= idx_j) THEN ! pass if k=i or k=j
            r_ik(:) = pos_i - pos_all(:,idx_k)
            CALL pbc(r_ik, box, hbox)
            d_ik = SQRT(SUM(r_ik**2))
            ! angle (k,i,j)
            N_ba = N_ba + 1
            dotprod = SUM(r_ij*r_ik)
            prod = d_ij*d_ik
            costheta_kij = dotprod/prod
            ! enforce cos(theta) >= -1
            IF (costheta_kij <= 0.0) THEN
              costheta_kij = DMAX1(-1.0_8,costheta_kij)
            END IF
            ! enforce cos(theta) <= 1
            IF (costheta_kij > 0.0) THEN
              costheta_kij = DMIN1(1.0_8,costheta_kij)
            END IF
            tetra = tetra + ABS(costheta_kij - costheta_tetra)
          END IF
        END DO
      END IF
    END DO
    tetra = tetra / N_ba
  END SUBROUTINE tetrahedrality
  
  
  !!!!!!!!!! TETRAHEDRALITY (ALL) !!!!!!!!!!
  SUBROUTINE tetrahedrality_all(idx_0, pos_0, pos_all, neigh, neigh_number, &
                                box, tetra)
    ! Parameters
    INTEGER(8), INTENT(in)  :: idx_0(:), neigh(:,:), neigh_number(:)
    REAL(8), INTENT(in)     :: pos_0(:,:), pos_all(:,:), box(:)
    REAL(8), INTENT(out)    :: tetra(SIZE(idx_0))
    ! Variables
    INTEGER(8) :: i, j, k, idx_i, idx_j, idx_k, nn_i, N_ba_i
    REAL(8)    :: r_ij(SIZE(box)), r_ik(SIZE(box)), d_ij, d_ik, hbox(SIZE(box))
    REAL(8)    :: dotprod, prod, costheta_kij, costheta_tetra
    ! Computation
    hbox = box / 2.0
    costheta_tetra = -0.333806859233771 ! cos(109.5°)
    tetra = 0.0
    DO i=1,SIZE(idx_0)
      idx_i = idx_0(i) + 1 ! python index shift
      nn_i = neigh_number(i)
      N_ba_i = 0
      ! first neighbor: j
      DO j=1,nn_i
        idx_j = neigh(i,j) + 1 ! python index shift
        IF (idx_j /= idx_i) THEN ! pass if j=i
          r_ij(:) = pos_0(:,i) - pos_all(:,idx_j)
          CALL pbc(r_ij, box, hbox)
          d_ij = SQRT(SUM(r_ij**2))
          ! second neighbor: k
          DO k=1,nn_i
            idx_k = neigh(i,k) + 1 ! python index shift
            IF (idx_k /= idx_i .AND. idx_k /= idx_j) THEN ! pass if k=i or k=j
              r_ik(:) = pos_0(:,i) - pos_all(:,idx_k)
              CALL pbc(r_ik, box, hbox)
              d_ik = SQRT(SUM(r_ik**2))
              ! angle (k,i,j)
              N_ba_i = N_ba_i + 1
              dotprod = SUM(r_ij*r_ik)
              prod = d_ij*d_ik
              costheta_kij = dotprod/prod
              ! enforce cos(theta) >= -1
              IF (costheta_kij <= 0.0) THEN
                costheta_kij = DMAX1(-1.0_8,costheta_kij)
              END IF
              ! enforce cos(theta) <= 1
              IF (costheta_kij > 0.0) THEN
                costheta_kij = DMIN1(1.0_8,costheta_kij)
              END IF
              tetra(i) = tetra(i) + ABS(costheta_kij - costheta_tetra)
            END IF
          END DO
        END IF
      END DO
      tetra(i) = tetra(i) / N_ba_i
    END DO
  END SUBROUTINE tetrahedrality_all


  !!!!!!!!!! PBC !!!!!!!!!!
  ! can probably be optimized further
  SUBROUTINE pbc_(r, box)
    REAL(8), INTENT(inout) :: r(:,:)
    REAL(8), INTENT(in)    :: box(:)
    INTEGER(8)             :: i
    REAL(8)                :: hbox(SIZE(box))
    hbox = box / 2.0    
    DO i=1,SIZE(box)
      WHERE (ABS(r(i,:)) > hbox(i))
        r(i,:) = r(i,:) - SIGN(box(i), r(i,:))
      END WHERE
    END DO
  END SUBROUTINE
  

  !!!!!!!!!! CARTESIAN TO SPHERICAL !!!!!!!!!!
  FUNCTION cartesian_to_spherical(r_xyz) RESULT(r_sph)
    REAL(8), INTENT(in) :: r_xyz(:,:)
    REAL(8)             :: xy(SIZE(r_xyz,2))
    REAL(8)             :: r_sph(SIZE(r_xyz,1),SIZE(r_xyz,2))
    xy = r_xyz(1,:)**2 + r_xyz(2,:)**2
    r_sph(1,:) = SQRT( xy + r_xyz(3,:)**2 ) ! module
    r_sph(2,:) = ATAN2( r_xyz(2,:), r_xyz(1,:) ) ! longitude
    r_sph(3,:) = ATAN2( SQRT(xy), r_xyz(3,:) ) ! latitude 
  END FUNCTION cartesian_to_spherical
  
  
  !!!!!!!! LEGENDRE FUNCTION !!!!!!!!!!
  FUNCTION plm(l, m, x) RESULT(plmx)
    INTEGER(8), INTENT(in) :: l ! degree
    INTEGER(8)             :: m ! order
    REAL(8), INTENT(in)    :: x(:) ! argument (must have |x| <= 1)
    REAL(8)                :: plmx(SIZE(x)) ! value of the associated Legendre function
    INTEGER(8)             :: i, ll
    REAL(8)                :: fact, pll(SIZE(x)), pmm(SIZE(x)), pmmp1(SIZE(x)), somx2(SIZE(x))
    LOGICAL                :: neg
    pll = 0.0_8
    neg = (m < 0)
    m = ABS(m)
    pmm = 1.0_8 ! compute P^m_m
    IF (m > 0) THEN
       somx2 = SQRT( (1.0_8 - x) * (1.0_8 + x) )
       fact = 1.0_8
       DO i=1,m
          pmm = -pmm * fact * somx2
          fact = fact + 2.0_8
       END DO
    END IF
    IF (l == m) THEN
      plmx = pmm
    ELSE
      pmmp1 = x * (2*m+1) * pmm ! compute P^m_m+1
      IF (l == m+1) THEN
        plmx = pmmp1
      ELSE ! compute P^m_l, l > m+1
        DO ll=m+2,l
          pll = (x * (2*ll-1) * pmmp1 - (ll+m-1) * pmm) / (ll-m)
          pmm = pmmp1
          pmmp1 = pll
        END DO
        plmx = pll
      END IF  
   END IF
   ! handle negative m
    IF (neg) THEN
      plmx = (-1.0_8)**m * REAL(factorial(l-m+1), 8) / REAL(factorial(l+m+1), 8) * plmx
    END IF
  END FUNCTION plm
  
  
  !!!!!!!!!! SPHERICAL HARMONICS !!!!!!!!!!
  FUNCTION ylm(l, m, theta, phi)
    INTEGER(8), INTENT(in) :: l, m
    REAL(8), INTENT(in)    :: theta(:), phi(:)
    COMPLEX(8)             :: ylm(SIZE(phi))
    REAL(8)                :: up, down
    INTEGER(16)            :: minus, plus
    minus = l - m
    plus = l + m
    up = (2*l+1) * REAL(factorial(minus+1), 8)
    down = 4.0*pi * REAL(factorial(plus+1), 8)
    ylm = SQRT(up/down) * plm(l, m, COS(phi)) * EXP( CMPLX(0.0, 1.0, KIND=8) * CMPLX(m*theta, 0.0, KIND=8) )
  END FUNCTION ylm
  

  !!!!!!!!!! COMPLEX VECTORS (INDIVIDUAL) !!!!!!!!!!
  FUNCTION qlm(l, neigh_i, pos_i, pos_all, box)
    INTEGER(8), INTENT(in) :: l, neigh_i(:)
    REAL(8), INTENT(in)    :: pos_i(:), pos_all(:,:), box(:)
    COMPLEX(8)             :: qlm(2*l+1), harm(SIZE(neigh_i))
    REAL(8)                :: r_xyz(3, SIZE(neigh_i)), r_sph(3, SIZE(neigh_i))
    INTEGER(8)             :: j, ni, m
    qlm(:) = (0.0, 0.0)
    ! r_ij (cartesian)
    DO j=1,SIZE(neigh_i)
      ni = neigh_i(j) + 1 ! python index shift 
      r_xyz(:,j) = pos_all(:,ni)
    END DO
    r_xyz(1,:) = r_xyz(1,:) - pos_i(1)
    r_xyz(2,:) = r_xyz(2,:) - pos_i(2)
    r_xyz(3,:) = r_xyz(3,:) - pos_i(3)
    CALL pbc_(r_xyz, box)
    ! r_ij (spherical)
    r_sph = cartesian_to_spherical(r_xyz)
    DO m=0,2*l
      harm = ylm(l, m-l, r_sph(2,:), r_sph(3,:))
      qlm(m+1) = qlm(m+1) + SUM(harm)
    END DO
    qlm = qlm / SIZE(neigh_i)
  END FUNCTION qlm

  !!!!!!!!!! COMPLEX VECTORS (ALL) !!!!!!!!!!
  SUBROUTINE qlm_all(l, neigh, neigh_number, pos_0, pos_all, box, q_lm)
    INTEGER(8), INTENT(in) :: l, neigh(:,:), neigh_number(:)
    REAL(8), INTENT(in)    :: pos_0(:,:), pos_all(:,:), box(:)
    COMPLEX(8)             :: q_lm(2*l+1, SIZE(neigh,1)), harm(SIZE(neigh,2))
    REAL(8)                :: r_xyz(3, SIZE(neigh,2)), r_sph(3, SIZE(neigh,2))
    INTEGER(8)             :: i, j, m, nn_i, k, idx_j
    DO i = 1,SIZE(pos_0,2)
       q_lm(:,i) = (0.0, 0.0)
       nn_i = neigh_number(i)
       DO j=1,nn_i
          ! TODO: here it would be better (j,i)
          idx_j = neigh(i,j) + 1 ! python index shift 
          r_xyz(:,j) = pos_all(:,idx_j)          
       END DO
       DO k = 1,SIZE(r_xyz,1)
          r_xyz(k,1:nn_i) = r_xyz(k,1:nn_i) - pos_0(k,i)
       END DO
       CALL pbc_(r_xyz(:,1:nn_i), box)
       ! r_ij (spherical)
       r_sph(:,1:nn_i) = cartesian_to_spherical(r_xyz(:,1:nn_i))
       DO m=0,2*l
          harm(1:nn_i) = ylm(l, m-l, r_sph(2,1:nn_i), r_sph(3,1:nn_i))
          q_lm(m+1,i) = q_lm(m+1,i) + SUM(harm(1:nn_i))
       END DO
       q_lm(:,i) = q_lm(:,i) / nn_i
    END DO
  END SUBROUTINE qlm_all


  !!!!!!!!!! ROTATIONAL INVARIANT OF ORDER l !!!!!!!!!!
  ! difference at the ~8th digit compared to Python 
  FUNCTION rotational_invariant(l, q_lm) RESULT(q_l)
    INTEGER(8), INTENT(in) :: l
    COMPLEX(8), INTENT(in) :: q_lm(:)
    REAL(8)                :: s, q_l
    s = REAL( SUM(q_lm * CONJG(q_lm)) )
    q_l = SQRT( 4.0*pi / REAL(2*l+1) * s )
  END FUNCTION rotational_invariant
    

  !!!!!!!!!! STEINHARDT (INDIVIDUAL) !!!!!!!!!!
  FUNCTION ql(l, neigh_i, pos_i, pos_all, box) RESULT(q_l)
    INTEGER(8), INTENT(in) :: l, neigh_i(:)
    REAL(8), INTENT(in)    :: pos_i(:), pos_all(:,:), box(:)    
    COMPLEX(8)             :: q_lm(2*l+1)
    REAL(8)                :: q_l
    q_lm = qlm(l, neigh_i, pos_i, pos_all, box)
    q_l = rotational_invariant(l, q_lm)  
  END FUNCTION ql
  

  !!!!!!!!!! STEINHARDT (ALL) !!!!!!!!!!
  SUBROUTINE ql_all(l, neigh, neigh_number, pos_0, pos_all, box, q_l)
    INTEGER(8), INTENT(in) :: l, neigh(:,:), neigh_number(:)
    REAL(8), INTENT(in)    :: pos_0(:,:), pos_all(:,:), box(:)    
    COMPLEX(8)             :: q_lm(2*l+1,SIZE(pos_0,2))
    REAL(8), INTENT(out)   :: q_l(SIZE(pos_0,2))
    INTEGER(8)             :: i
    CALL qlm_all(l, neigh, neigh_number, pos_0, pos_all, box, q_lm)
    DO i = 1,SIZE(pos_0,2)
       q_l(i) = rotational_invariant(l, q_lm(:,i))
    END DO
  END SUBROUTINE ql_all


!  !!!!!!!!!! COMPLEX VECTORS !!!!!!!!!!
!  ! For the generalization in Fortran of Lechner-Dellago 
!  FUNCTION qlm(l, neigh_i, pos_i, pos_j, box)
!    INTEGER(8), INTENT(in) :: l, neigh_i(:)
!    REAL(8), INTENT(in)    :: pos_i(:), pos_j(:,:), box(:)
!    COMPLEX(8)             :: qlm(2*l+1), harm
!    REAL(8)                :: r_xyz(3, SIZE(neigh_i)), r_sph(3, SIZE(neigh_i))
!    INTEGER(8)             :: j, n_neigh_i, ni, m
!    qlm(:) = (0.0, 0.0)
!    ! number of neighbors of i
!    n_neigh_i = 0
!    DO j=1,SIZE(neigh_i)
!      IF (neigh_i(j) /= -1) n_neigh_i = n_neigh_i + 1
!    END DO
!    ! r_ij (cartesian)
!    DO j=1,n_neigh_i
!      ni = neigh_i(j)+1 ! python index shift 
!      r_xyz(:,j) = pos_j(:,ni)
!    END DO
!    r_xyz(1,:) = r_xyz(1,:) - pos_i(1)
!    r_xyz(2,:) = r_xyz(2,:) - pos_i(2)
!    r_xyz(3,:) = r_xyz(3,:) - pos_i(3)
!    CALL pbc_(r_xyz, box)
!    ! r_ij (spherical)
!    r_sph = cartesian_to_spherical(r_xyz)
!    DO m=0,2*l
!      DO j=1,n_neigh_i
!        harm = ylm(l, m-l, r_sph(2,j), r_sph(3,j))
!        qlm(m+1) = qlm(m+1) + harm
!      END DO
!    END DO
!    qlm = qlm / n_neigh_i
!  END FUNCTION qlm


  !!!!!!!!!! AVERAGE COMPLEX VECTORS !!!!!!!!!!
!  FUNCTION qbarlm(l, neigh_i, neigh_neigh_i, pos_i, pos_j, box)
!    INTEGER(8), INTENT(in) :: l, neigh_i(:), neigh_neigh_i(:,:)
!    REAL(8), INTENT(in)    :: pos_i(:), pos_j(:,:), box(:)
!    INTEGER(8)             :: nbar_b, kn, k
!    COMPLEX(8)             :: q_lm_i(2*l+1), q_lm_k(SIZE(neigh_i), 2*l+1), qbarlm(2*l+1)
!    nbar_b = SIZE(neigh_i) + 1
!    q_lm_i = qlm(l, neigh_i, pos_i, pos_j, box)
!    DO kn=1,SIZE(neigh_i)
!      k = neigh_i(kn)+1 ! python index shift
!      q_lm_k(kn,:) = qlm(l, neigh_neigh_i(:,kn), pos_j(:,k), pos_j, box)
!    END DO
!    qbarlm = q_lm_i + SUM(q_lm_k, 1)
!    qbarlm = qbarlm / nbar_b
!  END FUNCTION qbarlm


  !!!!!!!!!! LECHNER-DELLAGO !!!!!!!!!!
!  FUNCTION qbarl(l, neigh_i, neigh_neigh_i, pos_i, pos_j, box) RESULT(qbar_l)
!    INTEGER(8), INTENT(in) :: l, neigh_i(:), neigh_neigh_i(:,:)
!    REAL(8), INTENT(in)    :: pos_i(:), pos_j(:,:), box(:)
!    COMPLEX(8)             :: qbar_lm(2*l+1)
!    REAL(8)                :: qbar_l
!    qbar_lm = qbarlm(l, neigh_i, neigh_neigh_i, pos_i, pos_j, box)
!    qbar_l = rotational_invariant(l, qbar_lm)
!  END FUNCTION qbarl
  

  !!!!!!!!!! SMOOTHED COMPLEX VECTORS (INDIVIDUAL) !!!!!!!!!!
  FUNCTION smoothed_qlm(l, neigh_i, pos_i, pos_all, spe_i, spe_all, cutoffs, pow, box) RESULT(q_lm)
    ! parameters
    INTEGER(8), INTENT(in) :: l, neigh_i(:), spe_i, spe_all(:), pow
    REAL(8), INTENT(in)    :: pos_i(:), pos_all(:,:), cutoffs(:,:), box(:)
    ! variables
    COMPLEX(8)             :: q_lm(2*l+1), harm(SIZE(neigh_i))
    REAL(8)                :: r_xyz(3, SIZE(neigh_i)), r_sph(3, SIZE(neigh_i))
    REAL(8)                :: d_ij(SIZE(neigh_i)), rc_ij, w_i(SIZE(neigh_i))
    INTEGER(8)             :: j, idx_j, m
    q_lm(:) = (0.0, 0.0)
    ! r_ij (cartesian)
    DO j=1,SIZE(neigh_i)
      idx_j = neigh_i(j) + 1 ! python index shift 
      r_xyz(:,j) = pos_all(:,idx_j)
    END DO
    r_xyz(1,:) = r_xyz(1,:) - pos_i(1)
    r_xyz(2,:) = r_xyz(2,:) - pos_i(2)
    r_xyz(3,:) = r_xyz(3,:) - pos_i(3)
    CALL pbc_(r_xyz, box)
    ! weights
    d_ij = SQRT(SUM(r_xyz**2, 1))
    DO j=1,SIZE(neigh_i)
      idx_j = neigh_i(j) + 1 ! python index shift
      rc_ij = cutoffs(spe_i, spe_all(idx_j))
      w_i(j) = EXP(-(d_ij(j) / rc_ij)**pow)
    END DO
    ! r_ij (spherical)
    r_sph = cartesian_to_spherical(r_xyz)
    DO m=0,2*l
      harm = ylm(l, m-l, r_sph(2,:), r_sph(3,:))
      q_lm(m+1) = q_lm(m+1) + DOT_PRODUCT(w_i, harm)
    END DO
    q_lm = q_lm / SUM(w_i)
  END FUNCTION smoothed_qlm


  !!!!!!!!!! SMOOTHED COMPLEX VECTORS (ALL) !!!!!!!!!!
  SUBROUTINE smoothed_qlm_all(l, neigh, neigh_number, pos_0, pos_all, spe, spe_all, cutoffs, pow, box, q_lm)
    ! TODO: fortran-wise it would be better to have neigh as (nmax, ndim) array
    ! parameters
    INTEGER(8), INTENT(in) :: l, neigh(:,:), neigh_number(:), spe(:), spe_all(:), pow
    REAL(8), INTENT(in)    :: pos_0(:,:), pos_all(:,:), cutoffs(:,:), box(:)
    ! variables
    COMPLEX(8), INTENT(inout) :: q_lm(2*l+1, SIZE(neigh,1))
    COMPLEX(8)                :: harm(SIZE(neigh,2))
    REAL(8)                   :: r_xyz(SIZE(pos_0,1), SIZE(neigh,2)), r_sph(SIZE(r_xyz,1), SIZE(r_xyz,2))
    ! All these arrays are cut down to size(neigh,2) (nmax)
    REAL(8)    :: d_ij(SIZE(neigh,2)), rc_ij, w_i(SIZE(neigh,2))
    INTEGER(8) :: i, j, k, idx_j, m, nn_i
    DO i = 1,SIZE(pos_0,2)
       q_lm(:,i) = (0.0, 0.0)
       nn_i = neigh_number(i)
       ! r_ij (cartesian)
       DO j=1,nn_i
          ! TODO: here it would be better (j,i)
          idx_j = neigh(i,j) + 1 ! python index shift 
          r_xyz(:,j) = pos_all(:,idx_j)          
       END DO
       DO k = 1,SIZE(r_xyz,1)
          r_xyz(k,1:nn_i) = r_xyz(k,1:nn_i) - pos_0(k,i)
       END DO
       CALL pbc_(r_xyz(:,1:nn_i), box)
       ! weights
       d_ij(1:nn_i) = SQRT(SUM(r_xyz(:, 1:nn_i)**2, 1))
       DO j = 1,nn_i
          idx_j = neigh(i,j) + 1 ! python index shift
          rc_ij = cutoffs(spe(i), spe_all(idx_j))
          w_i(j) = EXP(-(d_ij(j) / rc_ij)**pow)
       END DO
       ! r_ij (spherical)
       r_sph(:,1:nn_i) = cartesian_to_spherical(r_xyz(:,1:nn_i))
       DO m = 0,2*l
          harm(1:nn_i) = ylm(l, m-l, r_sph(2,1:nn_i), r_sph(3,1:nn_i))
          q_lm(m+1,i) = q_lm(m+1,i) + DOT_PRODUCT(w_i(1:nn_i), harm(1:nn_i))
       END DO
       q_lm(:,i) = q_lm(:,i) / SUM(w_i(1:nn_i))
    END DO
  END SUBROUTINE smoothed_qlm_all
  

  !!!!!!!!!! SMOOTHED STEINHARDT (INDIVIDUAL) !!!!!!!!!!
  FUNCTION smoothed_ql(l, neigh_i, pos_i, pos_all, spe_i, spe_all, box, cutoffs, pow) RESULT(q_l)
    INTEGER(8), INTENT(in) :: l, neigh_i(:), spe_i, spe_all(:), pow
    REAL(8), INTENT(in)    :: pos_i(:), pos_all(:,:), cutoffs(:,:), box(:)
    COMPLEX(8)             :: q_lm(2*l+1)
    REAL(8)                :: q_l
    q_lm = smoothed_qlm(l, neigh_i, pos_i, pos_all, spe_i, spe_all, cutoffs, pow, box)
    q_l = rotational_invariant(l, q_lm)  
  END FUNCTION smoothed_ql


  !!!!!!!!!! SMOOTHED STEINHARDT (ALL) !!!!!!!!!!
  SUBROUTINE smoothed_ql_all(l, neigh, neigh_number, pos_0, pos_all, spe, spe_all, box, cutoffs, pow, q_l)
    ! neigh: (neigh_max, npart)
    INTEGER(8), INTENT(in) :: l, neigh(:,:), neigh_number(:), spe(:), spe_all(:), pow
    REAL(8), INTENT(in)    :: pos_0(:,:), pos_all(:,:), cutoffs(:,:), box(:)
    COMPLEX(8)             :: q_lm(2*l+1,SIZE(pos_0,2))
    REAL(8), INTENT(out)   :: q_l(SIZE(pos_0,2))
    INTEGER                :: i
    CALL smoothed_qlm_all(l, neigh, neigh_number, pos_0, pos_all, spe, spe_all, cutoffs, pow, box, q_lm)
    DO i = 1,SIZE(pos_0,2)
       q_l(i) = rotational_invariant(l, q_lm(:,i))
    END DO
  END SUBROUTINE smoothed_ql_all


  !!!!!!!!!! DISTANCE-DEPENDENT COMPLEX VECTORS (INDIVIDUAL) !!!!!!!!!!
  FUNCTION radial_qlm(l, r, delta, exponent, neigh_i, pos_i, pos_all, box) RESULT(q_lmrd)
    ! parameters
    INTEGER(8), INTENT(in) :: l, neigh_i(:), exponent
    REAL(8), INTENT(in)    :: r, delta, pos_i(:), pos_all(:,:), box(:)
    ! variables
    COMPLEX(8)             :: q_lmrd(2*l+1), harm(SIZE(neigh_i))
    REAL(8)                :: r_xyz(3, SIZE(neigh_i)), r_sph(3, SIZE(neigh_i))
    REAL(8)                :: d_ij(SIZE(neigh_i)), Z(SIZE(neigh_i))
    INTEGER(8)             :: j, ni, m
    q_lmrd(:) = (0.0, 0.0)
    ! r_ij (cartesian)
    DO j=1,SIZE(neigh_i)
      ni = neigh_i(j) + 1 ! python index shift 
      r_xyz(:,j) = pos_all(:,ni)
    END DO
    r_xyz(1,:) = r_xyz(1,:) - pos_i(1)
    r_xyz(2,:) = r_xyz(2,:) - pos_i(2)
    r_xyz(3,:) = r_xyz(3,:) - pos_i(3)
    CALL pbc_(r_xyz, box)
    ! weights
    d_ij = SQRT(SUM(r_xyz**2, 1))
    Z = EXP(- 0.5d0 * ((d_ij - r) / delta)**exponent)
    ! r_ij (spherical)
    r_sph = cartesian_to_spherical(r_xyz)
    DO m=0,2*l
      harm = ylm(l, m-l, r_sph(2,:), r_sph(3,:))
      q_lmrd(m+1) = q_lmrd(m+1) + DOT_PRODUCT(Z, harm)
    END DO
    q_lmrd = q_lmrd / SUM(Z)
  END FUNCTION radial_qlm


  !!!!!!!!!! DISTANCE-DEPENDENT COMPLEX VECTORS (ALL) !!!!!!!!!!
  SUBROUTINE radial_qlm_all(l, r, delta, exponent, neigh, neigh_number, pos_0, &
                            pos_all, box, q_lmrd)
    ! parameters
    INTEGER(8), INTENT(in)    :: l, neigh(:,:), neigh_number(:), exponent
    REAL(8), INTENT(in)       :: r, delta, pos_0(:,:), pos_all(:,:), box(:)
    COMPLEX(8), INTENT(inout) :: q_lmrd(2*l+1,SIZE(neigh,1))
    ! variables
    COMPLEX(8)             :: harm(SIZE(neigh,2))
    REAL(8)                :: r_xyz(3, SIZE(neigh,2)), r_sph(3, SIZE(neigh,2))
    REAL(8)                :: d_ij(SIZE(neigh,2)), Z(SIZE(neigh,2))
    INTEGER(8)             :: i, j, nn_i, ni, dim, m
    DO i=1,SIZE(pos_0,2)
      q_lmrd(:,i) = (0.0, 0.0)
      nn_i = neigh_number(i)
      ! r_ij (cartesian)
      DO j=1,nn_i
        ni = neigh(i,j) + 1 ! python index shift
        r_xyz(:,j) = pos_all(:,ni)
      END DO
      DO dim=1,SIZE(r_xyz,1)
        r_xyz(dim,1:nn_i) = r_xyz(dim,1:nn_i) - pos_0(dim,i)
      END DO
      CALL pbc_(r_xyz(:,1:nn_i), box)
      ! weights
      d_ij(1:nn_i) = SQRT(SUM(r_xyz(:, 1:nn_i)**2, 1))
      Z = EXP(- 0.5d0 * ((d_ij - r) / delta)**exponent)
      ! r_ij (spherical)
      r_sph(:,1:nn_i) = cartesian_to_spherical(r_xyz(:,1:nn_i))
      DO m=0,2*l
        harm(1:nn_i) = ylm(l, m-l, r_sph(2,1:nn_i), r_sph(3,1:nn_i))
        q_lmrd(m+1,i) = q_lmrd(m+1,i) + DOT_PRODUCT(Z(1:nn_i), harm(1:nn_i))
      END DO
      q_lmrd(:,i) = q_lmrd(:,i) / SUM(Z(1:nn_i))
    END DO
  END SUBROUTINE radial_qlm_all


  !!!!!!!!!! ROTATIONAL INVARIANT OF DISTANCE-DEPENDENT BOP (INDIVIDUAL) !!!!!!!!!!
  FUNCTION radial_ql(l, r, delta, exponent, neigh_i, pos_i, pos_all, box) RESULT(q_lrd)
    INTEGER(8), INTENT(in) :: l, neigh_i(:), exponent
    REAL(8), INTENT(in)    :: r, delta, pos_i(:), pos_all(:,:), box(:)
    COMPLEX(8)             :: q_lmrd(2*l+1)
    REAL(8)                :: q_lrd
    q_lmrd = radial_qlm(l, r, delta, exponent, neigh_i, pos_i, pos_all, box)
    q_lrd = rotational_invariant(l, q_lmrd)  
  END FUNCTION radial_ql


  !!!!!!!!!! ROTATIONAL INVARIANT OF DISTANCE-DEPENDENT BOP (ALL) !!!!!!!!!!
  SUBROUTINE radial_ql_all(l, r, delta, exponent, neigh, neigh_number, pos_0, &
                         pos_all, box, q_lrd)
    INTEGER(8), INTENT(in) :: l, neigh(:,:), neigh_number(:), exponent
    REAL(8), INTENT(in)    :: r, delta, pos_0(:,:), pos_all(:,:), box(:)
    REAL(8), INTENT(out)   :: q_lrd(SIZE(pos_0,2))
    COMPLEX(8)             :: q_lmrd(2*l+1,SIZE(pos_0,2))
    INTEGER                :: i
    CALL radial_qlm_all(l, r, delta, exponent, neigh, neigh_number, pos_0, pos_all, box, q_lmrd)
    DO i=1,SIZE(pos_0,2)
      q_lrd(i) = rotational_invariant(l, q_lmrd(:,i))
    END DO
  END SUBROUTINE radial_ql_all
  

  !!!!!!!!!! COMPACTNESS !!!!!!!!!!
  SUBROUTINE compactness(pos, tetrahedra, radii, box, big_omega)
    ! parameters
    INTEGER(8), INTENT(in) :: tetrahedra(:,:)
    REAL(8), INTENT(in)    :: pos(:,:), radii(:), box(:)
    REAL(8), INTENT(OUT)   :: big_omega
    ! variables
    INTEGER(8)             :: tetra_idx, j, j_idx, k, k_idx
    REAL(8)                :: delta, norm, d_perfect, d_actual
    REAL(8)                :: r_jk(SIZE(box)), hbox(SIZE(box)), omega(SIZE(tetrahedra,2))
    ! computation
    hbox = box / 2.0
    DO tetra_idx=1,SIZE(tetrahedra,2)
      delta = 0.0
      norm = 0.0
      DO j=1,4
        j_idx = tetrahedra(j, tetra_idx) + 1 ! python index shift
        DO k=1,4
          k_idx = tetrahedra(k, tetra_idx) + 1 ! python index shift
          IF (k_idx > j_idx) THEN
            d_perfect = radii(j_idx) + radii(k_idx)
            r_jk(:) = pos(:,j_idx) - pos(:,k_idx)
            CALL pbc(r_jk, box, hbox)
            d_actual = SQRT(SUM(r_jk**2))
            delta = delta + ABS(d_actual - d_perfect)
            norm = norm + d_perfect
          END IF
        END DO
      END DO
      omega(tetra_idx) = delta / norm
    END DO
    big_omega = SUM(omega) / SIZE(omega)
  END SUBROUTINE compactness

  
END MODULE compute

