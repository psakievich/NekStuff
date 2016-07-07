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

      parameter(nfldm=5)
      parameter(lSYM=lx1*ly1*lz1)


      real pts      (ldim,lSYM)
     $               , fieldout (nfldm,lSYM)
     $               , dist     (lSYM)
     $               , rst      (lSYM*ldim)


      integer rcode(lSYM),elid(lSYM),proc(lSYM)

      common /scrcg/  pm1 (lx1,ly1,lz1,lelv) ! mapped pressure
      common /outtmp/ wrk (lx1*ly1*lz1*lelt,nfldm)
      character*3    prefix

      logical iffind

      integer icalld,npoints,npts,elmNum
      integer iEnd,iEndTotal
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


      if(nflds.ne.nfldm.and.nid.eq.0)write(6,*)"Error nflds ",nflds,
     $   nfldm
      call intpts_setup(-1.0,inth_hpts) ! use default tolerance
      elmNum=1 !initialize element
      iEnd=0 !set flag for end of elements to zero
      iEndTotal=0
c
c     BEGIN ELEMENT BASED LOOP
c
      do while (elmNum.le.nelt.and.iEnd.eq.0)
      if(nid.eq.0)then
         write(6,*),"Elm num",elmNum-1,"of",nelt
      end if
        call load_element(pts,npts,npoints,elmNum)
      if(nid.eq.0)then
         write(6,*),"Elm num",elmNum-1,"of",nelt
      end if
      !call exitt
      
      ! interpolate
        call findpts(inth_hpts,rcode,1,
     &                 proc,1,
     &                 elid,1,
     &                 rst,ndim,
     &                 dist,1,
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
      ! evaluate input field at given points
      do ifld = 1,nflds
         call findpts_eval(inth_hpts,fieldout(ifld,1),nfldm,
     &                     rcode,1,
     &                     proc,1,
     &                     elid,1,
     &                     rst,ndim,npts,
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

      if(elmNum.lt.nelt) then
        elmNum=elmNum+1
        if(nid.eq.0)write(6,*)"ElmNum Inc",elmNum
      else
        iEnd=1
      end if

      call gop(iEnd,iEndTotal,'*  ',1)
      if(nid.eq.0)then
         write(6,*),"Elm num",elmNum-1,"of",nelt,"iEnd equals",iend
      end if
      end do
      call prepost_map(1)  ! maps back axisymm arrays

      if(nio.eq.0) write(6,*) 'done :: swap points based on symmetry'

      return
      end
c-----------------------------------------------------------------------
      subroutine load_element(pts,npts,npoints,elmNum)
c     npts=local count; npoints=total count

      include 'SIZE'
      include 'TOTAL'
      !include 'PARALLEL'

      parameter (lt2=2*lx1*ly1*lz1*lelt)
      common /scrns/ xyz(ldim,lt2)
      common /scruz/ mid(lt2)  ! Target proc id
      integer i
      real    pts(ldim,npts)

      if(nid.eq.0)write(6,*)"Load elm ElmNum Inc",elmNum
      !load pnts
      do i=1,lx1*ly1*lz1
         pts(1,i)=xm1(i,1,1,elmNum)
         pts(2,i)=ym1(i,1,1,elmNum)
         pts(3,i)=1.0-zm1(i,1,1,elmNum)
      end do


      return
      end

