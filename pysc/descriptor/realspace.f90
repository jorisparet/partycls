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
 
!   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   !     Fortran implementation of the SANN algorithm                        !
!   !     van Meel, Filion, Valeriani and Frenkel November (2011)             !
!   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   ! TO FINISH
!   SUBROUTINE sann()
!     ! npart = total number of particles in the system
!     INTEGER(8) :: npart
!     !  m = tentative number of neighbours
!     INTEGER(8) :: i,j,k,m
!     ! countneighbors = number of neighbours of particle i
!     INTEGER(8) :: countneighbors(1000)
!     ! neighbor = list of neighbours of particles i
!     INTEGER(8) :: neighbor(1000,100)
!     ! sortneighbor = sorted neighbours 
!     INTEGER(8) :: sortneighbor(1000,100)
!     ! selectedneighbors = list of selected neighbours    
!     INTEGER(8) ::  selectedneighbors(1000,100)
!     ! Nb = final number of neighbours of particle i
!     INTEGER(8) :: Nb(1000)
!     ! edge of the simulation box
!     REAL(8) :: box(:)
!     REAL(8) :: hbox(:)
!     ! distance = list of distances between each 
!     ! neighbour of particle i and particle i 
!     REAL(8) ::  distance(1000,100)
!     ! distancesorted = sorted distances
!     REAL(8) :: distancesorted(1000,100)       
!     ! R(m) as in Eq.3 in the manuscript
!     REAL(8) :: rm,rm1
!     ! x,y,z component of every particle i
!     REAL(8) :: x(1000),y(1000),z(1000)
!     ! distance between particle i and particle j
!     REAL(8) :: d_ij
!     ! cutoff distance to identify all potential neighbours
!     REAL(8) :: rcutoff
!     
!     ! *** STEP 1 ***
!     ! first we identify the particles within a cutoff radius rcutoff
!     DO i=1,npart
!       ! loop over all particles different from i
!       DO j =1,npart
!         IF (j /= i) THEN
!           ! compute x,y,z component of the distance between particle i and j
!           r_ij(:) =  pos_1(:,j) - pos_0(:,i)
!           ! applying periodic boundary conditions
!           CALL pbc(r_ij, box, hbox)
!           ! compute distance d_ij  between particle i and j
!           d_ij = SQRT(SUM(r_ij**2))
!           ! compute distance d_ij  between particle i and j
!           d_ij = SQRT(dx*dx+dy*dy+dz*dz)
!           ! identify neighbours that are within a cutoff (rcutoff)
!           IF (d_ij < rcutoff) THEN
!             ! j is a neighbour of i
!             countneighbors(i) = countneighbors(i) + 1
!             ! build a list of neighbours
!             neighbor(i,countneighbors(i))= j 
!             ! create a list with the distance between i and j 
!             distance(i,countneighbors(i))= d_ij
!           END IF
!         END IF
!       END DO
!     END DO
!     
!     ! *** STEP 2 ***
!     ! for every particle i sort all (countneighbors) 
!     ! neighbours (neighbor) according to their 
!     ! distances (distance) and create  a new list of 
!     ! particle i's (sortneighbor)
!     ! and a new sorted list of distances (distancesorted)
!     DO i=1,npart
!          CALL sort(i,countneighbors,distance,neighbor,sortneighbor,distancesorted)
!     END DO
! 
!     DO i=1,npart
!       ! *** STEP 3 *** 
!       ! start with 3 neighbours
!       m = 3
!       ! *** STEP 4 ***
!       ! compute R(m) as in Eq.3 
!       rm = 0
!       DO k=1,m
!         rm = rm + distancesorted(i,k)
!       END DO
!       rm = rm/(m-2)
!       ! compute r(m+1)
!       DO j = 1,countneighbors(i)      
!         rm1 = 0
!         DO k=1,m
!            rm1 = rm1 + distancesorted(i,k)
!         END DO
!         rm1 = rm1/(m-2)
!         ! *** STEP 5 ***  
!         ! if rm > rm1     
!         IF (rm >= rm1) THEN     
!           rm = rm1
!           ! increase m
!           m = m+1
!         ELSE
!           ! *** STEP 6 ***
!           ! if rm < rm1, m is the final number of neighbours
!           EXIT
!         END IF
!       END DO
!       ! the final number of neighbours is m = Nb(i) 
!       ! and the neighbours are  selectedneighbors
!       Nb(i) = m
!       DO j=1,Nb(i)
!         selectedneighbors(i,j) = sortneighbor(i,j)
!       END DO
!     END DO
!       
!     RETURN
!     
!   END SUBROUTINE
  
  
END MODULE compute

