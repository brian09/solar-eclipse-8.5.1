      SUBROUTINE LOOP(FATHER,MOTHER,PERM,IERROR,MAXPEO,NPEO)
C
C     THIS SUBROUTINE DETERMINES WHETHER A PEDIGREE CONTAINS ANY
C     DIRECTED CYCLES.  IF SO, THEN SOMEONE IS HIS OWN ANCESTOR.
C     ON RETURN, IERROR IS 0 IF THERE ARE NO DIRECTED CYCLES.
C     OTHERWISE, IERROR IS THE NUMBER OF SOMEONE WHO IS HIS OWN
C     ANCESTOR.  THE PERMUTATION VECTOR PERM RETURNS WITH A REORDERING
C     OF THE PEDIGREE SO PARENTS PRECEDE THEIR CHILDREN.  FOR
C     ALGORITHMIC DETAILS SEE: E. LAWLER(1976) "COMBINATORIAL
C     OPTIMIZATION AND MATROIDS" HOLT, RINEHART, AND WINSTON, NEW YORK,
C     PAGE 32.  IF I IS A FOUNDER, THEN FATHER(I)=MOTHER(I)=0.
C     OTHERWISE, FATHER(I) IS A NUMBER BETWEEN 1 AND NPEO INDICATING
C     THE FATHER OF I.  DITTO FOR MOTHER.
C
      INTEGER FATHER(MAXPEO),MOTHER(MAXPEO),PERM(MAXPEO)        
C       
C     PUT ALL FOUNDERS AT THE START OF THE PERMUTATION VECTOR.  FOR
C     EACH REMAINING PERSON, PERM CONTAINS BOTH A CURRENT PERMUTATION
C     LOCATION AND HOW MANY OF HIS PARENTS HAVE BEEN ELIMINATED AS
C     POSSIBLE CANDIDATES FOR A DIRECTED CYCLE.
C
      IERROR=0
      M=1
      N=NPEO
      DO 10 I=1,NPEO
      IF (FATHER(I).EQ.0) THEN
      PERM(M)=I
      M=M+1
      ELSE
      PERM(N)=I+NPEO+NPEO
      N=N-1
      END IF
 10   CONTINUE
C
C     CHECK WHETHER ANYONE HAS BOTH PARENTS ELIMINATED.  SUCH A PERSON
C     CANNOT BELONG TO A DIRECTED CYCLE.  ELIMINATE THIS PERSON AND
C     COMPUTE HIS FINAL PERMUTATION LOCATION.  IF NO SUCH PERSON EXISTS,
C     THEN THE REMAINING PEOPLE ALL BELONG TO DIRECTED CYCLES.
C
      M=1
 50   DO 20 K=M,NPEO
 20   IF (PERM(K).LE.NPEO) GO TO 30
      IERROR=MOD(PERM(M)-1,NPEO)+1
      RETURN
 30   ISAVE=PERM(K)
      PERM(K)=PERM(M)
      PERM(M)=ISAVE
C
C     IF THIS IS THE LAST PERSON WE ARE DONE.
C
      M=M+1
      IF (M.GT.NPEO) RETURN
C
C     FIND THE CHILDRED OF THE ELIMINATED PERSON AND REDUCE THEIR
C     CURRENT PARENT COUNTS BY ONE.
C
      DO 40 I=M,NPEO
      IP=MOD(PERM(I)-1,NPEO)+1
      IF (FATHER(IP).EQ.ISAVE) PERM(I)=PERM(I)-NPEO
 40   IF (MOTHER(IP).EQ.ISAVE) PERM(I)=PERM(I)-NPEO
      GO TO 50
      END
