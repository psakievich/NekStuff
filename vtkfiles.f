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
      subroutine vtk_Interp
      !This subroutine outputs results from intpnts to a structured
      ! vtk file with parallel output (i.e. 1 file per processor)
      include 'SIZE'
      include 'TOTAL'

      parameter(nfldm=ldim+ldimt+1) !number of fields being interpolated

      !common block for findpts to put data into
      common /c_hptsr/ pts      (ldim,lhis)
     $               , fieldout (nfldm,lhis)
     $               , dist     (lhis)
     $               , rst      (lhis,ldim)

      common /c_hptsi/ rcode(lhis),elid(lhis),proc(lhis)
      common /scrcg/ pm1(lx1,ly1,lz1,lelv) !mapped pressure
      common /outtmp/ wrk(lx1*ly1*lz1*lelt,nfldm)

      !define local variables
      logical iffind
      integer icalled, npoints,npts
      save icalld,npoints,npts
      data icalld /0/
      data npoints /0/

      save inth_hpts !Not sure where this is defined, appears to be a
                     !handle based on intpts_setup

      nxyz = nx1*ny1*nz1
      ntot = nxyz*nelt
      nbuff = lhis

      if(nio.eq.0) write(6,*) 'dump vtk file'
      !0) CHECK IF ROUTINE NEEDS TO BE SETUP
      if(icalld.eq.0)then
         npts =lhis !number of points per processor

      !1) DETERMINE LOCAL RESPONSIBILITY POINTS |
      !2) CHECK THAT ENOUGH MEMORY IS AVAILABLE |
      !3) ASSIGN POINTS                         |
      !4) SETUP INTERPOLATE ROUTINES            |
      !                                        \|/
         call vtk_hpts_in(pts,npts,npoints) !populate pts and npoints
         call intpts_setup(-1.0,inth_hpts) !use default tolerance
      endif
      !5) MAP PRESSURE
      !6) PACK WORKING ARRAY
      !7) INTERPOLATE
      !8) CHECK RETURN CODES
      !9) EVALUTAE FIELD AT GIVEN POINTS
      !10) WRITE VTK_PARALLEL FILE ON RANK0
      !11) WRITE VTK FILE(S) ON ALL RANKS
      !12) MAP PRESSURE BACK
      end subroutine
C ---------------------------------------------------------------------
c-----------------------------------------------------------------------
      subroutine vtk_hpts_in(pts,npts,npoints)
c                        npts=local count; npoints=total count

      include 'SIZE'
      include 'PARALLEL'

      parameter (lt2=2*lx1*ly1*lz1*lelt)
      common /scrns/ xyz(ldim,lt2)
      common /scruz/ mid(lt2)  ! Target proc id
      common/hpts_to_elm/NELGH,NXH,NYH,NZH !sizes from hpts_fld
      real    pts(ldim,npts)

c    I think that if this conditional is false the routine
c    puts all of the points on processor 0
      if (lt2.gt.npts) then

         call ps_buffer_in(xyz,npp,npoints,lt2) !lt2 is the size of buffer
         if(npoints.gt.np*npts) then
           if(nid.eq.0)write(6,*)'ABORT in hpts(): npoints > NP*lhis!!'
           if(nid.eq.0)write(6,*)'Change SIZE: ',np,npts,npoints
           call exitt
         endif
         if(npoints.ne.NELGH*NXH*NYH*NZH)then
           if(nid.eq.0)write(6,*)'HPNTS dosnt match the given dims'
     $                     ,npoints,NELGH*NXH*NYH*NZH
           call exitt
         endif

         npmax = (npoints/npts)
         if(mod(npoints,npts).eq.0) npmax=npmax+1

         if(nid.gt.0.and.npp.gt.0) then
          npts_b = lt2*(nid-1)               ! # pts  offset(w/o 0)
          nprc_b = npts_b/npts               ! # proc offset(w/o 0)

          istart = mod(npts_b,npts)          ! istart-->npts pts left
          ip     = nprc_b + 1                ! PID offset
          icount = istart                    ! point offset
         elseif(nid.eq.0) then
          npts0   = mod1(npoints,lt2)        ! Node 0 pts
          npts_b  = npoints - npts0          ! # pts before Node 0
          nprc_b  = npts_b/npts

          istart  = mod(npts_b,npts)
          ip      = nprc_b + 1
          icount  = istart
         endif

         do i =1,npp
            icount = icount + 1
            if(ip.gt.npmax) ip = 0
            mid(i) = ip
            if (icount.eq.npts) then
               ip     = ip+1
               icount = 0
            endif
         enddo

         call crystal_tuple_transfer
     &      (cr_h,npp,lt2,mid,1,pts,0,xyz,ldim,1)

         call copy(pts,xyz,ldim*npp)
      else
         call ps_buffer_in(pts,npp,npoints,npts)
      endif
      npts = npp


      return
      end
c-----------------------------------------------------------------------
