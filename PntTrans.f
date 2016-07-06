c-----------------------------------------------------------------------
c SUBROUTINES DEVELOPED BY PHIL SAKIEVICH
c FOR TRANSFORMATIONS OF GRID
c   
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
      subroutine SymFlip()
c    ********************************************
c    ***** MODIFIED VERSION OF HPTS IN REPO******
c     flips all fields about the mid plane of the simulation
c    ********************************************
C
c     evaluate velocity, temperature, pressure and ps-scalars 
c     for list of points (read from hpts.in) and dump results
c     into a file (hpts.out).
c     note: read/write on rank0 only 
c
c     ASSUMING LHIS IS MAX NUMBER OF POINTS TO READ IN ON ONE PROCESSOR

      include 'SIZE'
      include 'TOTAL'

      parameter(nfldm=ldim+ldimt+1)
      parameter(lSYM=lx1*ly1*lz1)


      common /SYMc_hptsr/ SYMPts      (ldim,lSYM)
     $               , SYMfieldout (nfldm,lSYM)
     $               , SYMdist     (lSYM)
     $               , SYMrst      (lSYM*ldim)


      common /SYMc_hptsi/ SYMrcode(lSYM),SYMelid(lSYM),SYMproc(lSYM)

      common /scrcg/  pm1 (lx1,ly1,lz1,lelv) ! mapped pressure
      common /outtmp/ wrk (lx1*ly1*lz1*lelt,nfldm)
      character*3    prefix

      logical iffind

      integer icalld,npoints,npts,elmNum
      save    icalld,npoints,npts
      data    icalld  /0/
      data    npoints /0/

      save    inth_hpts

      nxyz  = nx1*ny1*nz1
      ntot  = nxyz*nelt 
      nbuff = lhis      ! point to be read in on 1 proc.
      npts = nxyz
      if(nio.eq.0) write(6,*) 'swap points based on symmetry'
      call prepost_map(0)  ! maps axisymm and pressure
      ! pack working array
      ! modified to dump out corrdinates as well
      nflds = ndim
      if(ifvo) then
        call copy(wrk(1,1),vx,ntot)
        call copy(wrk(1,2),vy,ntot)
        if(if3d) call copy(wrk(1,3),vz,ntot)
        nflds = ndim
      endif
      if(ifpo) then
        nflds = nflds + 1
        call copy(wrk(1,nflds),pm1,ntot)
      endif
      if(ifto) then
        nflds = nflds + 1
        call copy(wrk(1,nflds),t,ntot)
      endif
      do i = 1,ldimt
         if(ifpsco(i)) then
           nflds = nflds + 1
           call copy(wrk(1,nflds),T(1,1,1,1,i+1),ntot)
         endif
      enddo

c     BEGIN ELEMENT BASED LOOP
      elmNum=1 !initialize element
      do while (elmNum.le.nelt)

        call load_element(pts,npts,npoints,elmNum)
        call intpts_setup(-1.0,inth_hpts) ! use default tolerance

      
      ! interpolate
      if(icalld.eq.0) then
        call findpts(inth_hpts,SYMrcode,1,
     &                 SYMproc,1,
     &                 SYMelid,1,
     &                 SYMrst,ndim,
     &                 SYMdist,1,
     &                 pts(1,1),ndim,
     &                 pts(2,1),ndim,
     &                 pts(3,1),ndim,npts)
      
        do i=1,npts
           ! check return code 
           if(rcode(i).eq.1) then
             if (dist(i).gt.1e-12) then
                nfail = nfail + 1
                IF (NFAIL.LE.5) WRITE(6,'(a,1p4e15.7)') 
     &     ' WARNING: point on boundary or outside the mesh xy[z]d^2:'
     &     ,(pts(k,i),k=1,ndim),dist(i)
             endif   
           elseif(rcode(i).eq.2) then
             nfail = nfail + 1
             if (nfail.le.5) write(6,'(a,1p3e15.7)') 
     &        ' WARNING: point not within mesh xy[z]: !',
     &        (pts(k,i),k=1,ndim)
           endif
        enddo
        icalld = 1
      endif
      if(nflds.ne.nfldm.and.nid.eq.0)write(6,*)"Error nflds ",nflds,
     $   nfldm
      ! evaluate input field at given points
      do ifld = 1,nflds
         call findpts_eval(inth_hpts,fieldout(ifld,1),nfldm,
     &                     SYMrcode,1,
     &                     SYMproc,1,
     &                     SYMelid,1,
     &                     SYMrst,ndim,npts,
     &                     wrk(1,ifld))
      enddo

      
      !Write results back to the current element
      do i=1,nxyz
         vx(i,1,1,elmNum)=-fieldout(1,i)
         vy(i,1,1,elmNum)=-fieldout(2,i)
         vz(i,1,1,elmNum)=-fieldout(3,i)
         pm1(i,1,1,elmNum)=-fieldout(4,i)
         t(i,1,1,elmNum,1)=1.0-fieldout(5,i)
      end do
      end do
      call prepost_map(1)  ! maps back axisymm arrays

      if(nio.eq.0) write(6,*) 'done :: swap points based on symmetry'

      return
      end
c-----------------------------------------------------------------------
      subroutine load_element(pts,npts,npoints,elmNum)
c     npts=local count; npoints=total count

      include 'SIZE'
      include 'PARALLEL'

      parameter (lt2=2*lx1*ly1*lz1*lelt)
      common /scrns/ xyz(ldim,lt2)
      common /scruz/ mid(lt2)  ! Target proc id
      integer elmNum,i
      real    pts(ldim,npts)

      !load pnts
      do i=1,lx1*ly1*lz1
         pts(1,i)=xm1(i,1,1,elmNum)
         pts(2,i)=ym1(i,1,1,elmNum)
         pts(3,i)=1.0-zm1(i,1,1,elmNum)
      end do


      return
      end
