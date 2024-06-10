      program Brownian
c
c this is to illustrate Brownian dynamics in 1D
c
      integer nstep
      real*8 fact1,fact2,kbt,dt,t
      real*8 rgauss
      external rgauss
      real*8 delx
      real*8 k
      real  rho(10000)

      iseed = 837388130
      t = 0.0
      dt = 0.005
      dx = 0.01
      k =  1.0
      xwindow = -2.5
      x = xwindow 
      diffusion = 1.0
      kBT = 0.5915
      nstep = 1000000
      xmin = -5.0
      xmax = 5.0
      do i=1,10000
      rho(i) = 0.0
      enddo
      do istep=1,nstep

c        call ENERGY(fx,x)
c energy = 0.1*x**4-x**2
         Fx = -0.1*4*x**3+2*x         
c        Fx = Fx -k*(x-xwindow)

         fact1 = diffusion*dt/kBT
         delx = Fx*fact1
         fact2 = sqrt(2*diffusion*dt)
         x = x + delx + fact2*rgauss(iseed)
         t = t + dt

         if((x.gt.xmin).and.(x.lt.xmax))then
         ix = int((x-xmin)/dx) + 1
         rho(ix) = rho(ix) + 1
         endif
         write(10,'(2f16.3)') t,x

      enddo
      tot = 0.0
      do i = 1,10000
      tot = tot + rho(i)
c      write(*,*) 'tot',tot
      enddo

      do i=1,10000
      x = xmin+(i-1)*dx+dx/2
      if((x.gt.xmin).and.(x.lt.xmax))then
      write(11,'(f16.3,f16.8)') x ,rho(i)/tot 
      endif
      enddo

      END


c-------------------------------------------------------------------------------

      FUNCTION RGAUSS(ISEED)
c     random gaussian numbers
      real*8     twopi
      parameter (twopi=2*3.141592658979D0)
      real*8 rgauss, random
      integer iseed
      external random
      rgauss = cos(random(iseed)*twopi)*SQRT(-2*LOG(random(iseed)))
      write(*,*) random(iseed),random(iseed)
      return
      end


      FUNCTION RANDOM(ISEED)
      REAL*8 RANDOM
C     RANDOM NUMBER GENERATOR: UNIFORM DISTRIBUTION (0,1)
C     ISEED: SEED FOR GENERATOR. ON THE FIRST CALL THIS HAS TO
C     HAVE A VALUE IN THE EXCLUSIVE RANGE (1, 2147483647)
C     AND WILL BE REPLACED BY A NEW VALUE TO BE USED IN
C     FOLLOWING CALL.
C
C     REF: Lewis, P.A.W., Goodman, A.S. & Miller, J.M. (1969)
C     "Pseudo-random number generator for the System/360", IBM
C     Systems Journal 8, 136.
C
C     This is a "high-quality" machine independent generator.
C     INTEGERS are supposed to be 32 bits or more.
C     The same algorithm is used as the basic IMSL generator.
C
C     Author: Lennart Nilsson
C
      INTEGER ISEED
      REAL*8 DSEED,DIVIS,DENOM,MULTIP
      DATA  DIVIS/2147483647.D0/
      DATA  DENOM /2147483711.D0/
      DATA  MULTIP/16807.D0/
C
      IF(ISEED.LE.1) ISEED=314159
      DSEED=MULTIP*ISEED
      DSEED=MOD(DSEED,DIVIS)
      RANDOM=DSEED/DENOM
      ISEED=DSEED
C
      RETURN
      END
