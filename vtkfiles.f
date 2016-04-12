C ---------------------------------------------------------------------
C     SUBROUTINES DEVELOPED BY PHIL SAKIEVICH
C     EMAIL: PSAKIEVI@ASU.EDU
C
C     WRITE INTERPOLATED DATA TO VTK_STRUCTURED_GRID XML FORMAT
C     IN PARALLEL
C ---------------------------------------------------------------------
      subroutine vtk_Test
      !Test that the file compiles correctly
      !with nek make routines
      write(6,*),"VTK ROUTINES ARE COMPILING"
      end subroutine vtk_Test
C ---------------------------------------------------------------------
      subroutine vtk_PartitionCylinder(nR,nTheta,nZ)
      !This subroutine partitions the domain and sets up the intpnts
      !on each processor.  The r-z plane is divided up between
      !processors, but all of theta is maintained on each processor.
      !i.e. each processor contains "rings" of data
      include 'SIZE'
      include 'TOTAL'
      integer nR, nTheta, nZ
      end subroutine
C ---------------------------------------------------------------------
