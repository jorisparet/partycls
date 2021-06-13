MODULE compute

  IMPLICIT NONE

  ! Constant Pi
  REAL(8),  PARAMETER :: pi  = 4.0_8 * ATAN(1.0_8)
  
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
          bin = FLOOR(SQRT(d_ij_sq) / dr) + FLOOR(rmin / dr) + 1
          IF (bin >= 1 .AND. bin <= nbins) hist(i, bin) = hist(i, bin) + 1
        END IF
      END DO
    END DO
  END SUBROUTINE

  
  SUBROUTINE angular_histogram(idx_i, pos_i, pos_1, neigh_i, box, nbins, dtheta, hist)
    ! Parameters
    INTEGER(8), INTENT(in)  :: idx_i, neigh_i(:), nbins
    REAL(8), INTENT(in)     :: pos_i(:), pos_1(:,:), box(:), dtheta
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
      IF (idx_j /= idx_i) THEN ! pass if j=i
        r_ij(:) = pos_i - pos_1(:,idx_j)
        CALL pbc(r_ij, box, hbox)
        d_ij = SQRT(SUM(r_ij**2))
        ! Second neighbor: k
        DO k=1,n_neigh_i
          idx_k = neigh_i(k) + 1 ! python index shift
          IF (idx_k /= idx_i .AND. idx_k /= idx_j) THEN ! pass if k=i or k=j
            r_ik(:) = pos_i - pos_1(:,idx_k)
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
  END SUBROUTINE


  !!!!!!!!!! PBC !!!!!!!!!!
  ! Can be optimized
  SUBROUTINE pbc_(r, box)
    REAL(8), INTENT(inout) :: r(:,:)
    REAL(8), INTENT(in)    :: box(:)
    INTEGER(8)             :: i, j
    REAL(8)                :: hbox(SIZE(box))
    hbox = box/2.0    
    DO i=1,SIZE(box)
      DO j=1,SIZE(r,2)
        IF (ABS(r(i,j)) > hbox(i)) THEN
          r(i,j) = r(i,j) - SIGN(box(i), r(i,j)) 
        END IF    
      END DO    
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

  !!!!!!!!!! FACTORIAL !!!!!!!!!!
  RECURSIVE FUNCTION factorial(n) RESULT(fact)
    INTEGER(8) :: n, fact
    IF (n == 0) THEN
       fact = 1
    ELSE
       fact = n * factorial(n-1)
    END IF
  END FUNCTION factorial
  
  !!!!!!!! LEGENDRE FUNCTION !!!!!!!!!!
  FUNCTION plm(l, m, x) RESULT(plmx)
    INTEGER(8), INTENT(in) :: l ! degree
    INTEGER(8)             :: m ! order
    REAL(8), INTENT(in)    :: x ! argument (must have |x| <= 1)
    REAL(8)                :: plmx ! value of the associated Legendre function
    INTEGER(8)             :: i, ll
    REAL(8)                :: fact, pll, pmm, pmmp1, somx2
    LOGICAL                :: neg
    neg = (m < 0)
    m = ABS(m)
    pmm = 1.0 ! compute P^m_m
    IF (m > 0) THEN
       somx2 = SQRT( (1.0-x)*(1.0+x) )
       fact = 1.0
       DO i=1,m
          pmm = -pmm * fact * somx2
          fact = fact + 2.0
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
      plmx = (-1.0)**m * REAL(factorial(l-m))/REAL(factorial(l+m)) * plmx
    END IF
  END FUNCTION plm

  !!!!!!!!!! SPHERICAL HARMONICS !!!!!!!!!!
  FUNCTION ylm(l, m, theta, phi)
    INTEGER(8), INTENT(in) :: l, m
    REAL(8), INTENT(in)    :: theta, phi
    COMPLEX(8)             :: ylm
    REAL(8)                :: up, down
    up = (2*l+1)*factorial(l-m)
    down = 4.0*pi*factorial(l+m)
    !ylm = SQRT(up/down) * plm(l, m, COS(theta)) * EXP(CMPLX(0.0, 1.0)*CMPLX(REAL(m*phi), 0.0))
    ylm = SQRT(up/down) * plm(l, m, COS(phi)) * EXP(CMPLX(0.0, 1.0)*CMPLX(REAL(m*theta), 0.0))
  END FUNCTION ylm
  
  !!!!!!!!!! COMPLEX VECTORS !!!!!!!!!!
  FUNCTION qlm(l, neigh_i, pos_i, pos_j, box)
    INTEGER(8), INTENT(in) :: l, neigh_i(:)
    REAL(8), INTENT(in)    :: pos_i(:), pos_j(:,:), box(:)
    COMPLEX(8)             :: qlm(2*l+1), harm
    REAL(8)                :: r_xyz(3, SIZE(neigh_i)), r_sph(3, SIZE(neigh_i))
    INTEGER(8)             :: j, ni, m
    qlm(:) = (0.0, 0.0)
    ! r_ij (cartesian)
    DO j=1,SIZE(neigh_i)
      ni = neigh_i(j)+1 ! python index shift 
      r_xyz(:,j) = pos_j(:,ni)
    END DO
    r_xyz(1,:) = r_xyz(1,:) - pos_i(1)
    r_xyz(2,:) = r_xyz(2,:) - pos_i(2)
    r_xyz(3,:) = r_xyz(3,:) - pos_i(3)
    CALL pbc_(r_xyz, box)
    ! r_ij (spherical)
    r_sph = cartesian_to_spherical(r_xyz)
    DO m=0,2*l
      DO j=1,SIZE(neigh_i)
        harm = ylm(l, m-l, r_sph(2,j), r_sph(3,j))
        qlm(m+1) = qlm(m+1) + harm
      END DO
    END DO
    qlm = qlm / SIZE(neigh_i)
  END FUNCTION qlm

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
  
  !!!!!!!!!! ROTATIONAL INVARIANT OF ORDER l !!!!!!!!!!
  ! difference at the ~8th digit compared to Python 
  FUNCTION rotational_invariant(l, q_lm) RESULT(q_l)
    INTEGER(8), INTENT(in) :: l
    COMPLEX(8), INTENT(in) :: q_lm(:)
    REAL(8)                :: s, q_l
    s = REAL( SUM(q_lm * CONJG(q_lm)) )
    q_l = SQRT( 4.0*pi / REAL(2*l+1) * s )
  END FUNCTION rotational_invariant
    
  !!!!!!!!!! STEINHARDT !!!!!!!!!!
  FUNCTION ql(l, neigh_i, pos_i, pos_j, box) RESULT(q_l)
    INTEGER(8), INTENT(in) :: l, neigh_i(:)
    REAL(8), INTENT(in)    :: pos_i(:), pos_j(:,:), box(:)    
    COMPLEX(8)             :: q_lm(2*l+1)
    REAL(8)                :: q_l
    q_lm = qlm(l, neigh_i, pos_i, pos_j, box)
    q_l = rotational_invariant(l, q_lm)  
  END FUNCTION ql

  !!!!!!!!!! AVERAGE COMPLEX VECTORS !!!!!!!!!!
  FUNCTION qbarlm(l, neigh_i, neigh_neigh_i, pos_i, pos_j, box)
    INTEGER(8), INTENT(in) :: l, neigh_i(:), neigh_neigh_i(:,:)
    REAL(8), INTENT(in)    :: pos_i(:), pos_j(:,:), box(:)
    INTEGER(8)             :: nbar_b, kn, k
    COMPLEX(8)             :: q_lm_i(2*l+1), q_lm_k(SIZE(neigh_i), 2*l+1), qbarlm(2*l+1)
    nbar_b = SIZE(neigh_i) + 1
    q_lm_i = qlm(l, neigh_i, pos_i, pos_j, box)
    DO kn=1,SIZE(neigh_i)
      k = neigh_i(kn)+1 ! python index shift
      q_lm_k(kn,:) = qlm(l, neigh_neigh_i(:,kn), pos_j(:,k), pos_j, box)
    END DO
    qbarlm = q_lm_i + SUM(q_lm_k, 1)
    qbarlm = qbarlm / nbar_b
  END FUNCTION qbarlm

  !!!!!!!!!! LECHNER-DELLAGO !!!!!!!!!!
  FUNCTION qbarl(l, neigh_i, neigh_neigh_i, pos_i, pos_j, box) RESULT(qbar_l)
    INTEGER(8), INTENT(in) :: l, neigh_i(:), neigh_neigh_i(:,:)
    REAL(8), INTENT(in)    :: pos_i(:), pos_j(:,:), box(:)
    COMPLEX(8)             :: qbar_lm(2*l+1)
    REAL(8)                :: qbar_l
    qbar_lm = qbarlm(l, neigh_i, neigh_neigh_i, pos_i, pos_j, box)
    qbar_l = rotational_invariant(l, qbar_lm)
  END FUNCTION qbarl
  
    
  SUBROUTINE nearest_neighbors(idx_i, idx_1, pos_i, pos_1, spe_i, spe_1, pairs, box, cutoffs, neigh_i)
    ! Parameters
    INTEGER(8), INTENT(in)     :: idx_i, idx_1(:)
    REAL(8), INTENT(in)        :: pos_i(:), pos_1(:,:)
    INTEGER(8), INTENT(in)     :: spe_i, spe_1(:)
    INTEGER(8), INTENT(in)     :: pairs(:,:)
    REAL(8), INTENT(in)        :: box(:)
    REAL(8), INTENT(in)        :: cutoffs(:)
    INTEGER(8), INTENT(out)    :: neigh_i(100) ! max. number of neighbors is assumed to be 100
    ! Variables
    INTEGER(8)                 :: np1, j, idxj, n_neigh_i, spe_j, line
    REAL(8)                    :: r_ij(SIZE(box)), dij_sq, rcut, rcut_sq, hbox(SIZE(box))
    ! Computation
    hbox = box / 2.0
    np1 = SIZE(spe_1)
    neigh_i = -1 ! set all neighbors to -1
    n_neigh_i = 0
    DO j=1,np1
      idxj = idx_1(j)
      IF (idxj /= idx_i) THEN
        spe_j = spe_1(j)
        ! find appropriate cutoff for pair (i,j)
        line = 1
        DO WHILE (pairs(line,1) /= spe_i .OR. pairs(line,2) /= spe_j)
          line = line + 1
        END DO
        rcut = cutoffs(line)
        rcut_sq = rcut*rcut
        ! test j as a neighbor of i
        r_ij(:) = pos_i - pos_1(:,j)
        CALL pbc(r_ij, box, hbox)
        dij_sq = SUM(r_ij**2)
        IF (dij_sq <= rcut_sq) THEN
          n_neigh_i = n_neigh_i + 1
          neigh_i(n_neigh_i) = j-1 ! python index shift
        END IF
      END IF
    END DO
  END SUBROUTINE nearest_neighbors
  

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !     Fortran implementation of the SANN algorithm                        !
  !     van Meel, Filion, Valeriani and Frenkel November (2011)             !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE sann(pos_i, pos_all, idx_i, idx_all, idx_1, rcutoff, box, selectedneighbors_1)
    ! Parameters
    REAL(8), INTENT(in)     :: pos_i(:), pos_all(:,:) ! position of i and all other particles
    INTEGER(8), INTENT(in)  :: idx_i, idx_all(:), idx_1(:) ! indices of i and all other particles
    REAL(8), INTENT(in)     :: rcutoff ! cutoff distance to identify all potential neighbours 
    REAL(8), INTENT(in)     :: box(:) ! edge of the simulation box  
    INTEGER(8), INTENT(out) :: selectedneighbors_1(250) ! list of selected neighbours in group = 1
    ! Variables
    INTEGER(8) :: j, k, m ! m = tentative number of neighbours
    INTEGER(8) :: np ! total number of particles, particles in group=1
    INTEGER(8) :: countneighbors ! countneighbors = number of neighbours of particle i
    INTEGER(8) :: neighbor(250) ! list of neighbours of particles i
    INTEGER(8) :: sortneighbor(250) ! sorted 
    INTEGER(8) :: selectedneighbors(250) ! list of selected neighbours
    INTEGER(8) :: Nb ! final number of neighbours of particle i
    INTEGER(8) :: idx_j ! index of particle j
    REAL(8)    :: distance(250) ! distance = list of distances between each 
    REAL(8)    :: distancesorted(250) ! distancesorted = sorted distances       
    REAL(8)    :: rm, rm1 ! R(m) as in Eq.3 in the manuscript
    REAL(8)    :: r_ij(SIZE(box)) ! vector between i and j
    REAL(8)    :: d_ij ! distance between particle i and particle j
    ! Computation
    np = SIZE(idx_all)
    countneighbors = 0
    selectedneighbors_1 = -1
    ! *** STEP 1 ***
    ! first we identify the particles within a cutoff radius rcutoff
    ! loop over all particles different from i
    DO j =1,np
      idx_j = idx_all(j)
      IF (idx_j /= idx_i) THEN
        ! compute x,y,z component of the distance between particle i and j
        r_ij(:) =  pos_all(:,j) - pos_i
        ! applying periodic boundary conditions
        CALL pbc(r_ij, box, box/2.0)
        ! compute distance d_ij  between particle i and j
        d_ij = SQRT(SUM(r_ij**2))
        ! identify neighbours that are within a cutoff (rcutoff)
        IF (d_ij < rcutoff) THEN
          ! j is a neighbour of i
          countneighbors = countneighbors + 1
          ! build a list of neighbours
          neighbor(countneighbors) = j 
          ! create a list with the distance between i and j 
          distance(countneighbors) = d_ij
        END IF
      END IF
    END DO
    
    ! *** STEP 2 ***
    ! sort all (countneighbors)
    ! neighbours (neighbor) according to their 
    ! distances (distance) and create  a new list of 
    ! particle i's (sortneighbor)
    ! and a new sorted list of distances (distancesorted)
    CALL sort(countneighbors, distance, neighbor, sortneighbor, distancesorted)
    
    ! *** STEP 3 *** 
    ! start with 3 neighbours
    m = 3
    ! *** STEP 4 ***
    ! compute R(m) as in Eq.3 
    rm = 0
    DO k=1,m
      rm = rm + distancesorted(k)
    END DO
    rm = rm/(m-2)
    ! compute r(m+1)
    DO j = 1,countneighbors
      rm1 = 0
      DO k=1,m
         rm1 = rm1 + distancesorted(k)
      END DO
      rm1 = rm1/(m-2)
      ! *** STEP 5 ***  
      ! if rm > rm1     
      IF (rm >= rm1) THEN     
        rm = rm1
        ! increase m
        m = m+1
      ELSE
        ! *** STEP 6 ***
        ! if rm < rm1, m is the final number of neighbours
        EXIT
      END IF
    END DO
    ! the final number of neighbours is m = Nb(i) 
    ! and the neighbours are  selectedneighbors
    Nb = m
    DO j=1,Nb
      selectedneighbors(j) = sortneighbor(j) - 1 ! python index shift
    END DO
    ! neighbors that are not in group=1 must be removed from the list
    k = 1
    DO j=1,Nb
      IF (ANY(idx_1 == idx_all(selectedneighbors(j)))) THEN
        selectedneighbors_1(k) = selectedneighbors(j)
        k = k + 1
      END IF
    END DO
  END SUBROUTINE
 

  SUBROUTINE sort(countneighbors, distance, neighbor, sortneighbor, distancesorted)
    ! Parameters
    INTEGER(8), INTENT(in)    :: countneighbors
    INTEGER(8), INTENT(inout) :: neighbor(:)
    REAL(8), INTENT(inout)    :: distance(:)
    INTEGER(8), INTENT(inout) :: sortneighbor(:)
    REAL(8), INTENT(inout)    :: distancesorted(:)
    ! Variables
    INTEGER(8) :: i, imin, j, n_tmp
    REAL(8)    :: d_tmp
    ! Computation
    DO i=1,countneighbors
      imin = i
      DO j=i+1,countneighbors
        IF (distance(j) < distance(imin)) THEN
          imin = j
        END IF
        d_tmp = distance(i)
        n_tmp = neighbor(i)
        distance(i) = distance(imin)
        neighbor(i) = neighbor(imin)
        distance(imin) = d_tmp
        neighbor(imin) = n_tmp
        distancesorted(i) = distance(i)
        sortneighbor(i) = neighbor(i)
      END DO
    END DO
  END SUBROUTINE  
  
END MODULE compute

