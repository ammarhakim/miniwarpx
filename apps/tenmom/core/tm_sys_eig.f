c
c     ==================================================================
      subroutine tm_sys_eig(meqn,qli,qri,wr,vl,vr)
c     ==================================================================
c
c     Computes flux eigensystem for 10 moment equations
c     
      implicit none
c     INTENT(IN)
      integer meqn
      double precision qli(10),qri(10)
c     INTENT(OUT)
      double precision wr(10),vl(10,10),vr(10,10)

c     local variables/arrays
c     stuff needed for LAPACK routines
      integer i,j,info
c      integer iwork,ipiv(meqn)
      double precision q(10),dfdq(10,10)
      double precision wi(10),work(10*100)
      double precision vll(10,10),vlvr(10,10)
c     
      double precision r,u,v,w,wxx,wxy,wxz,wyy,wyz,wzz

c     iwork = meqn*100          ! LAPACK work space size
      
c     compute average state at interface
      do i = 1,meqn
         q(i) = 0.5d0*(qli(i)+qri(i))
      end do

c     compute primitive variables
      r = q(1)
      u = q(2)/q(1)
      v = q(3)/q(1)
      w = q(4)/q(1)
      wxx = (q(5)-r*u*u)/r
      wxy = (q(6)-r*u*v)/r
      wxz = (q(7)-r*u*w)/r
      wyy = (q(8)-r*v*v)/r
      wyz = (q(9)-r*v*w)/r
      wzz = (q(10)-r*w*w)/r

c     compute flux jacobian
      do i = 1,meqn
         do j = 1,meqn
c     ----- initialize everything to zero
            dfdq(i,j) = 0.d0
         end do
      end do

c     ++ row 1:
      dfdq(1,2) = 1.d0
c     ++ row 2
      dfdq(2,5) = 1.d0

c     ++ row 3
      dfdq(3,6) = 1.d0
c     ++ row 4
      dfdq(4,7) = 1.d0
c     ++ row 5
      dfdq(5,1) = u**3 - 3.d0*u*wxx
      dfdq(5,2) = -3.d0*u**2 + 3.d0*wxx
      dfdq(5,5) = 3.d0*u
c     ++ row 6
      dfdq(6,1) = v*u**2 - 2.d0*u*wxy - v*wxx
      dfdq(6,2) = -2.d0*u*v + 2.d0*wxy
      dfdq(6,3) = -u**2 + wxx
      dfdq(6,5) = v
      dfdq(6,6) = 2.d0*u
c     ++ row 7
      dfdq(7,1) = w*u**2 - 2.d0*u*wxz - w*wxx
      dfdq(7,2) = -2.d0*u*w + 2.d0*wxz
      dfdq(7,4) = -u**2 + wxx
      dfdq(7,5) = w
      dfdq(7,7) = 2.d0*u
c     ++ row 8
      dfdq(8,1) = u*v**2 - u*wyy - 2.d0*v*wxy
      dfdq(8,2) = -v**2 + wyy
      dfdq(8,3) = -2.d0*u*v + 2.d0*wxy
      dfdq(8,6) = 2.d0*v
      dfdq(8,8) = u
c     ++ row 9
      dfdq(9,1) = u*v*w - u*wyz - v*wxz - w*wxy
      dfdq(9,2) = -v*w + wyz
      dfdq(9,3) = -u*w + wxz
      dfdq(9,4) = -u*v + wxy
      dfdq(9,6) = w
      dfdq(9,7) = v
      dfdq(9,9) = u
c     ++ row 10
      dfdq(10,1) = u*w**2 - u*wzz - 2.d0*w*wxz
      dfdq(10,2) = -w**2 + wzz
      dfdq(10,4) = -2.d0*u*w + 2.d0*wxz
      dfdq(10,7) = 2.d0*w
      dfdq(10,10) = u

c     compute eigenvalues and right eigenvectors
      call dgeev('V','V',meqn,dfdq,meqn,wr,wi,vll,meqn,vr,meqn,work,
     $     100,info)

c     find the product of the left and right eigenvectors. This should
c     be a diagonal matrix
      call dgemm('T','N',meqn,meqn,meqn,1.d0,vll,meqn,
     $     vr,meqn,
     $     0.d0,
     $     vlvr,meqn)

c     renormalize left eigenvectors to give correct normalization
c     (vl'*vr = unit matrix)
      do i = 1,meqn
c         write(*,*) i, vlvr(i,i)
         do j = 1,meqn
            vll(i,j) = vll(i,j)/vlvr(j,j)
         end do
      end do

c     transpose left eigenvector matrix,
      do i = 1,meqn
         do j = 1,meqn
            vl(i,j) = vll(j,i)
         end do
      end do

      return
c     -- subroutine tm_sys_eig ------------------------------------------
      end
