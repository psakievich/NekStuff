C***********************************************************************
C     ROUTINES FOR PERFORMING FFT's INSIDE NEK5000
C     REQUIRES FFTW (www.fftw.org), BUILT AND TESTED WITH v3.3.4
C
C     CODE WRITEN BY PHIL SAKIEVICH (psakievi@asu.edu)
C
C     UTILIZES DEFINITIONS IN fftw3.f by Matteo Frigo and Steven Johnson
C      (http://people.sc.fsu.edu/~jburkardt/f77_src/fftw3/fftw3.html)
C
C     INCLUDE FILE:  fftw3.f
C     LINK LIBS:    -lfftw3 -lw
C***********************************************************************

C     OVERVIEW:
C     These subroutines are designed to allow one to set up a set of
C     points on a given processor get their values from somewhere in the
C     parallel envirnoment using intpts and then perform FFT's on them
C     locally using routines from fftw3.3.4.
C
C     It is the users responsibility to ensure that the:
C
C     1) The way the points are defined are compatible with the FFT they
C        intend to perform.

C     SUBROUTINE MYFFT
      subroutine MyFFT()
      include 'SIZE'
      include 'TOTAL'
      include 'MYFFT'


      integer nFFTSetup  !variable to determine if setup has been called
      data nFFTSetup /0/ !initialize value to zero
      save nFFTSetup     !save value between subsequent calls

      ! 1) Make sure a valid number of processors are present
      if(nFFTProcs<np) then
         if(nid.eq.0)write(6,*),"ERROR FFT: FFT Procs< Total processors"
         call exitt()
      end if

      ! 2) Set up points
      if(nFFTSetup.eq.0) then
           call FFTDefinePoints()
           !call FFTFindPoints()
      end if

      ! 3) Perform Interpolation
           !call FFTInterpPoints()
      ! 4) If desired write to file

      return
      end

C     SUBROUTINE DEFINE POINTS
C     Users should use this to define the sampling points they want for
C     their FFT's. This example will be for points in a cylinder with
C     FFT in the theta direction
      subroutine FFTDefinePoints()
      include 'SIZE'
      include 'TOTAL'
      include 'MYFFT'

      parameter,real::PI=4.0*atan(1.0)
      integer i,j,k,ii

      real dR=1.0/nFFTlx1, dTheta=2.0*PI/nFFTly1, dZ=1.0/nFFTlz1

      ! Good idea to zero out any thing that won't be using an FFT
      if(nid.gt.nFFTprocs) then
       do i=1,nFFTtotal
          rFFTpts(1,i)=0.0
          rFFTpts(2,i)=0.0
          if(if3d) rFFTpts(3,1)=0.0
       end do
      else !Initialize fields for processors of interest
       do i=1,nFFTlx1
          do j=1,nFFTly1
             do k=1,nFFTlz1
                ii=k+j*nFFTlz1+i*nFFTlz1*nFFTly1
                rFFTpts(1,ii)=(i-1)*dR
                rFFTpts(2,ii)=(j-1)*dTheta
                rFFTpts(3,ii)=(k-1)*dZ
             end do
          end do
       end do
      end if

      return
      end
C     SUBROUTINE FIND PNTS
C     SUBROUTINE CREATE FFT PLANS
C     SUBROUTINE INTERPOLATE POINTS
C     SUBROUTINE PERFORM FFT
C     SUBROUTINE PRINT FFT DATA TO TEXT
