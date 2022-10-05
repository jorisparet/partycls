MODULE nearest_neighbors

  IMPLICIT NONE

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

  
  !!!!!!!!!! FIXED-CUTOFFS ALL PARTICLES !!!!!!!!!!
  SUBROUTINE fixed_cutoffs_all(pos, spe, box, cutoffs_sq, nn, neigh)
    ! Parameters
    REAL(8), INTENT(in)       :: pos(:,:)
    INTEGER(8), INTENT(in)    :: spe(:)
    REAL(8), INTENT(in)       :: box(:)
    REAL(8), INTENT(in)       :: cutoffs_sq(:,:)
    INTEGER(8), INTENT(inout) :: nn(:), neigh(:,:)
    ! Variables
    INTEGER(8) :: i, j, spe_i, spe_j
    REAL(8)    :: hbox(SIZE(box)), r_ij(SIZE(box)), dij_sq
    ! Computation
    hbox = box / 2.0
    DO i=1,SIZE(pos,2)
      spe_i = spe(i)
      DO j=i+1,SIZE(pos,2)
        spe_j = spe(j)
        r_ij = pos(:,i) - pos(:,j)
        CALL pbc(r_ij, box, hbox)
        dij_sq = DOT_PRODUCT(r_ij, r_ij)
        IF (dij_sq < cutoffs_sq(spe_i, spe_j)) THEN
          nn(i) = nn(i) + 1
          nn(j) = nn(j) + 1
          neigh(i, nn(i)) = j - 1 ! python index shift
          neigh(j, nn(j)) = i - 1 ! python index shift
        END IF
      END DO
    END DO
  END SUBROUTINE fixed_cutoffs_all


  !!!!!!!!!! FIXED-CUTOFFS DISTINCT SETS !!!!!!!!!!
  SUBROUTINE fixed_cutoffs_distinct(idx_0, idx_1, pos_0, pos_1, spe_0, spe_1, box, cutoffs_sq, nn_0, neigh_0)
    ! Parameters
    INTEGER(8), INTENT(in)    :: idx_0(:), idx_1(:)
    REAL(8), INTENT(in)       :: pos_0(:,:), pos_1(:,:)
    INTEGER(8), INTENT(in)    :: spe_0(:), spe_1(:)
    REAL(8), INTENT(in)       :: box(:)
    REAL(8), INTENT(in)       :: cutoffs_sq(:,:)
    INTEGER(8), INTENT(inout) :: nn_0(:), neigh_0(:,:)
    ! Variables
    INTEGER(8) :: i, j, idx_i, idx_j, spe_i, spe_j
    REAL(8)    :: hbox(SIZE(box)), r_ij(SIZE(box)), dij_sq
    ! Computation
    hbox = box / 2.0
    DO i=1,SIZE(idx_0)
      idx_i = idx_0(i)
      spe_i = spe_0(i)
      DO j=1,SIZE(idx_1)
        idx_j = idx_1(j)
        spe_j = spe_1(j)
        IF (idx_j /= idx_i) THEN
          r_ij = pos_0(:,i) - pos_1(:,j)
          CALL pbc(r_ij, box, hbox)
          dij_sq = DOT_PRODUCT(r_ij, r_ij)
          IF (dij_sq < cutoffs_sq(spe_i, spe_j)) THEN
            nn_0(i) = nn_0(i) + 1
            neigh_0(i, nn_0(i)) = idx_j
          END IF
        END IF
      END DO
    END DO
  END SUBROUTINE fixed_cutoffs_distinct
  
  
  !!!!!!!!!! SANN !!!!!!!!!!
  ! van Meel, Filion, Valeriani and Frenkel November (2011)
  SUBROUTINE sann(pos_i, pos_all, idx_i, idx_all, rcutoff, box, selectedneighbors)
    ! Parameters
    REAL(8), INTENT(in)     :: pos_i(:), pos_all(:,:) ! position of i and all other particles
    INTEGER(8), INTENT(in)  :: idx_i, idx_all(:) ! indices of i and all other particles
    REAL(8), INTENT(in)     :: rcutoff ! cutoff distance to identify all potential neighbours 
    REAL(8), INTENT(in)     :: box(:) ! edge of the simulation box  
    INTEGER(8), INTENT(out) :: selectedneighbors(500) ! list of selected neighbours in group = 1
    ! Variables
    INTEGER(8) :: j, k, m ! m = tentative number of neighbours
    INTEGER(8) :: np ! total number of particles, particles in group=1
    INTEGER(8) :: countneighbors ! countneighbors = number of neighbours of particle i
    INTEGER(8) :: neighbor(500) ! list of neighbours of particles i
    INTEGER(8) :: sortneighbor(500) ! sorted 
    INTEGER(8) :: Nb ! final number of neighbours of particle i
    INTEGER(8) :: idx_j ! index of particle j
    REAL(8)    :: distance(500) ! distance = list of distances between each 
    REAL(8)    :: distancesorted(500) ! distancesorted = sorted distances       
    REAL(8)    :: rm, rm1 ! R(m) as in Eq.3 in the manuscript
    REAL(8)    :: r_ij(SIZE(box)) ! vector between i and j
    REAL(8)    :: d_ij ! distance between particle i and particle j
    ! Computation
    np = SIZE(idx_all)
    countneighbors = 0
    selectedneighbors = -1
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
  
END MODULE nearest_neighbors
