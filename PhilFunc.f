c-----------------------------------------------------------------------
c SUBROUTINES DEVELOPED BY PHIL SAKIEVICH
c FOR RAYLEIGH-BENARD CONVECTION ANALYSIS
c   
c-----------------------------------------------------------------------
c
      subroutine ps_Test
c     TEST that the file compiles correctly 
c     with the nek make routines
       write(6,*),"PHIL ROUTINES ARE WORKING"      
      end subroutine ps_Test 
c-----------------------------------------------------------------------
      subroutine ps_T_Diss(psi)
c     this subroutine computes the thermal dissipation grad T:grad T at each
c     grid point
c
c     
      include 'SIZE'
      include 'TOTAL'
c    
      integer psi !which passive scalar to use to store variable
      real dTx(lx1,ly1,lz1,lelt),
     $     dTy(lx1,ly1,lz1,lelt),
     $     dTz(lx1,ly1,lz1,lelt)
c
      nt=nx1*ny1*nz1*nelt
c     compute gradients
      call gradm1(dTx,dTy,dTz,T(1,1,1,1,1))
c
c     zero out the passive scalar
      call rzero(T(1,1,1,1,psi),nt)
c
c     
      call add2col2(T(1,1,1,1,psi),dTx,dTx,nt)
      call add2col2(T(1,1,1,1,psi),dTy,dTy,nt)
      call add2col2(T(1,1,1,1,psi),dTz,dTz,nt)
      return
      end subroutine ps_T_Diss
c-----------------------------------------------------------------------
      subroutine ps_KE_Diss(psi)
c     this subroutine computes (grad U+ grad U^T)^2 at each grid point
c     *note additional scaling will be required based on dimensional form
c     of equations and puts dissipation in passive scalar #psi
c
c
      include 'SIZE'
      include 'TOTAL'
      integer psi
c     derivatives of velocity field
      real ddux(lx1,ly1,lz1,lelt),
     $     dduy(lx1,ly1,lz1,lelt),
     $     dduz(lx1,ly1,lz1,lelt),
     $     ddvx(lx1,ly1,lz1,lelt),
     $     ddvy(lx1,ly1,lz1,lelt),
     $     ddvz(lx1,ly1,lz1,lelt),
     $     ddwx(lx1,ly1,lz1,lelt),
     $     ddwy(lx1,ly1,lz1,lelt),
     $     ddwz(lx1,ly1,lz1,lelt),
     $     work(lx1,ly1,lz1,lelt)
      integer nv
c
      nv=nx1*ny1*nz1*nelv
      nt=nx1*ny1*nz1*nelt
      call rzero(T(1,1,1,1,psi),nt) ! zero out epsilon
c   
c    compute velocity gradients
      call gradm1(ddux,dduy,dduz,vx)
      call gradm1(ddvx,ddvy,ddvz,vy)
      call gradm1(ddwx,ddwy,ddwz,vz)
c    sum up terms that contribute to dissipation
      call add2col2(T(1,1,1,1,psi),ddux,ddux,nt) !ux^2
      call add2col2(T(1,1,1,1,psi),ddvy,ddvy,nt) !vy^2
      call add2col2(T(1,1,1,1,psi),ddwz,ddwz,nt) !wz^2
      call add2col2(T(1,1,1,1,psi),dduy,ddvx,nt) 
      call add2col2(T(1,1,1,1,psi),ddvz,ddwy,nt) 
      call add2col2(T(1,1,1,1,psi),ddwx,dduz,nt) 
c   terms above contribute twice
      call cmult(T(1,1,1,1,psi),2.0,nt)
c   add the rest of the terms
      call add2col2(T(1,1,1,1,psi),dduy,dduy,nt)
      call add2col2(T(1,1,1,1,psi),dduz,dduz,nt)
      call add2col2(T(1,1,1,1,psi),ddvx,ddvx,nt)
      call add2col2(T(1,1,1,1,psi),ddvz,ddvz,nt)
      call add2col2(T(1,1,1,1,psi),ddwx,ddwx,nt)
      call add2col2(T(1,1,1,1,psi),ddwy,ddwy,nt)
c 
c      call rzero(T(1,1,1,1,psi),nt) ! zero out epsilon
      return        
      end subroutine ps_KE_Diss
c-----------------------------------------------------------------------
      subroutine ps_Dissipation(eps,n)
c     this subroutine computes (grad U+ grad U^T)^2 at each grid point
c     *note additional scaling will be required based on dimensional form
c     of equations
c
c     parameters: eps- real array for storing the dissipation
c                  n - integer for size of eps
c
      include 'SIZE'
      include 'TOTAL'
      integer n
      real eps(n)
c     derivatives of velocity field
      real dux(lx1,ly1,lz1,lelv),
     $     duy(lx1,ly1,lz1,lelv),
     $     duz(lx1,ly1,lz1,lelv),
     $     dvx(lx1,ly1,lz1,lelv),
     $     dvy(lx1,ly1,lz1,lelv),
     $     dvz(lx1,ly1,lz1,lelv),
     $     dwx(lx1,ly1,lz1,lelv),
     $     dwy(lx1,ly1,lz1,lelv),
     $     dwz(lx1,ly1,lz1,lelv),
     $     work(lx1,ly1,lz1,lelv)
      integer nv
c
      nv=nx1*ny1*nz1*nelv
      if(n.ne.nv)then
         write(6,*)"Error in divergence input size",n,nv
         return
      endif
      call rzero(eps,n) ! zero out epsilon
c   
c    compute velocity gradients
      call gradm1(dux,duy,duz,vx,nv)
      call gradm1(dvx,dvy,dvz,vy,nv)
      call gradm1(dwx,dwy,dwz,vz,nv)
c    sum up terms that contribute to dissipation
      call add2col2(eps,dux,dux,nv) !ux^2
      call add2col2(eps,dvy,dvy,nv) !vy^2
      call add2col2(eps,dwz,dwz,nv) !wz^2
      call add2col2(eps,duy,dvx,nv) 
      call add2col2(eps,dvz,dwy,nv) 
      call add2col2(eps,dwx,duz,nv) 
c   terms above contribute twice
      call cmult(eps,2.0,nv)
c   add the rest of the terms
      call add2col2(eps,duy,duy,nv)
      call add2col2(eps,duz,duz,nv)
      call add2col2(eps,dvx,dvx,nv)
      call add2col2(eps,dvz,dvz,nv)
      call add2col2(eps,dwx,dwx,nv)
      call add2col2(eps,dwy,dwy,nv)
c 
      return        
      end subroutine ps_Dissipation
c-----------------------------------------------------------------------
c
c
c-----------------------------------------------------------------------
      subroutine ps_hpts(prefix)
c    ********************************************
c    ***** MODIFIED VERSION OF HPTS IN REPO******
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

      parameter(nfldm=2*ldim+ldimt+1)


      common /c_hptsr/ pts      (ldim,lhis)
     $               , fieldout (nfldm,lhis)
     $               , dist     (lhis)
     $               , rst      (lhis*ldim)


      common /c_hptsi/ rcode(lhis),elid(lhis),proc(lhis)

      common /scrcg/  pm1 (lx1,ly1,lz1,lelv) ! mapped pressure
      common /outtmp/ wrk (lx1*ly1*lz1*lelt,nfldm)
      character*3    prefix

      logical iffind

      integer icalld,npoints,npts
      save    icalld,npoints,npts
      data    icalld  /0/
      data    npoints /0/

      save    inth_hpts

      nxyz  = nx1*ny1*nz1
      ntot  = nxyz*nelt 
      nbuff = lhis      ! point to be read in on 1 proc.

      if(nio.eq.0) write(6,*) 'dump history points'

      if(icalld.eq.0) then
        npts  = lhis      ! number of points per proc
        call ps_hpts_in(pts,npts,npoints) !npts and npoints are blank???
        call intpts_setup(-1.0,inth_hpts) ! use default tolerance
      endif


      call prepost_map(0)  ! maps axisymm and pressure

      ! pack working array
      ! modified to dump out corrdinates as well
      nflds = ndim
      if(ifvo) then
        call copy(wrk(1,ndim+1),vx,ntot)
        call copy(wrk(1,ndim+2),vy,ntot)
        if(if3d) call copy(wrk(1,ndim+3),vz,ntot)
        nflds = ndim+ndim
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
      
      ! interpolate
      if(icalld.eq.0) then
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
        icalld = 1
      endif
      if(nflds.ne.nfldm.and.nid.eq.0)write(6,*)"Error nflds ",nflds,
     $   nfldm
      ! evaluate input field at given points
      do ifld = ndim+1,nflds
         call findpts_eval(inth_hpts,fieldout(ifld,1),nfldm,
     &                     rcode,1,
     &                     proc,1,
     &                     elid,1,
     &                     rst,ndim,npts,
     &                     wrk(1,ifld))
      enddo
      !copy coordinates
      do i=1,ndim
        do ii=1,npts
          fieldout(i,ii)=pts(i,ii)
        enddo
      enddo
      
      ! write interpolation results to file
      call ps_hpts_out(fieldout,nflds,nfldm,npoints,nbuff)
c      call ps_hpts_out_fld(prefix,fieldout,nflds,nfldm,
c     $                                npoints,nbuff)

      call prepost_map(1)  ! maps back axisymm arrays

      if(nio.eq.0) write(6,*) 'done :: dump history points'

      return
      end
c-----------------------------------------------------------------------
      subroutine ps_buffer_in(buffer,npp,npoints,nbuf)
        
      include 'SIZE'
      include 'PARALLEL'

      common/hpts_to_elm/NELGH,NXH,NYH,NZH !sizes from hpts_fld
      real    buffer(ldim,nbuf)  

      ierr = 0
c    read in the total number of history points from hpts.in
      if(nid.eq.0) then
        write(6,*) 'reading ps_hpts.in'
        open(50,file='ps_hpts.in',status='old',err=100)
        read(50,*,err=100) npoints,NELGH,NXH,NYH,NZH
        goto 101
 100    ierr = 1
 101    continue
      endif
c   check all processors for an error
      ierr=iglsum(ierr,1)
      if(ierr.gt.0) then
        write(6,*) 'Cannot open ps_hpts.in in 
     $                    subroutine hpts()'
        call exitt
      endif
c    send total number of points to all processors and 
c    check to see if there is enough memory allocated      
      call bcast(npoints,isize)
      call bcast(nelgh,isize)
      call bcast(nxh,isize)
      call bcast(nyh,isize)
      call bcast(nzh,isize)
      if(npoints.gt.lhis*np) then
        if(nid.eq.0) write(6,*) 'ABORT: Too many pts to read in hpts()!'
        call exitt
      endif
      if(nid.eq.0) write(6,*) 'found ', npoints, ' points'

c    nbuf=2*lx1*ly1*lz1*lelt
      npass =  npoints/nbuf +1  !number of passes to cover all pts
      n0    =  mod(npoints,nbuf)!remainder 
      if(n0.eq.0) then
         npass = npass-1
         n0    = nbuf
      endif

      len = wdsize*ndim*nbuf
c    put all processors except processor 0 into receive mode
      if (nid.gt.0.and.nid.lt.npass) msg_id=irecv(nid,buffer,len)
      call nekgsync
c    read in data with processor 0 from hts and send to other
c    processors      
      npp=0  
      if(nid.eq.0) then
        i1 = nbuf
        do ipass = 1,npass
           if(ipass.eq.npass) i1 = n0
           do i = 1,i1
              read(50,*) (buffer(j,i),j=1,ndim) 
           enddo
           if(ipass.lt.npass)call csend(ipass,buffer,len,ipass,0)
        enddo
        close(50)
        npp = n0
c        open(50,file='hpts.out')!,status='new')
c        write(50,'(A)') 
c     &      '# time  vx  vy  [vz]  pr  T  PS1  PS2  ...'
      elseif (nid.lt.npass)  then !processors receiving data
        call msgwait(msg_id)
        npp=nbuf
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine ps_hpts_in(pts,npts,npoints) 
c                        npts=local count; npoints=total count

      include 'SIZE'
      include 'PARALLEL'

      parameter (lt2=2*lx1*ly1*lz1*lelt)
      common /scrns/ xyz(ldim,lt2)
      common /scruz/ mid(lt2)  ! Target proc id
      common/hpts_to_elm/NELGH,NXH,NYH,NZH !sizes from hpts_fld
      real    pts(ldim,npts)

c    I think that if this conditional is true the routine
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
      subroutine ps_hpts_out(fieldout,nflds,nfldm,npoints,nbuff)
c    ********************************************
c    *** MODIFIED VERSION OF HPTS_OUT IN REPO****
c    ********************************************

      include 'SIZE'
      include 'TOTAL'
      common/hpts_to_elm/NELGH,NXH,NYH,NZH !sizes from hpts_fld
      real buf(nfldm,nbuff),fieldout(nfldm,nbuff)
      character*80 filename
      character*1 excode(30)
      integer iFileNum
      save iFileNum
      data iFileNum /0/
      
      len = wdsize*nfldm*nbuff


      npass = npoints/nbuff + 1
      il = mod(npoints,nbuff)
      if(il.eq.0) then
         il = nbuff
         npass = npass-1
      endif
      !setup the header
      call BLANK(EXCODE,30)
         IF(IFXYO) then
            EXCODE(1)='X'
            EXCODE(2)=' '
            EXCODE(3)='Y'
            EXCODE(4)=' '
            i = 5
            IF(IF3D) THEN
              EXCODE(i)  ='Z'
              EXCODE(i+1)=' '
              i = i + 2
            ENDIF
         ENDIF
         IF(IFVO) then
            EXCODE(i)  ='U'
            EXCODE(i+1)=' '
            i = i + 2
         ENDIF
         IF(IFPO) THEN
           EXCODE(i)='P'
           EXCODE(i+1)=' '
           i = i + 2
         ENDIF
         IF(IFTO) THEN
           EXCODE(i)='T '
           EXCODE(i+1)=' '
           i = i + 1
         ENDIF
         do iip=1,ldimt1
            if (ifpsco(iip)) then
              write(excode(iip+I)  ,'(i1)') iip
              write(excode(iip+I+1),'(a1)') ' '
              i = i + 1
            endif
         enddo
 

        if(nid.eq.0)then
          write(filename,"('udfpnts.fld',I2.2)")iFileNum+1
          open(unit=50,file=filename)!,status='new')
          WRITE(50,'(4i4,1pe14.7,I5,1X,30A1,1X,A12)')
     $      NELGH,NXH,NYH,NZH,TIME,iFileNum,(EXCODE(I),I=1,30),
     $           'NELT,NX,NY,N'
          CDRROR=0.0
          WRITE(50,'(6G11.4)')(CDRROR,I=1,NELGH)   ! dummy 

        endif
        call nekgsync
      do ipass = 1,npass

        if(ipass.lt.npass) then
          if(nid.eq.0) then
            call crecv(ipass,buf,len)
            do ip = 1,nbuff
              write(50,'(1p20E14.6)'),
c     &         (pts(i,ip), i=1,ndim),
     &         (buf(i,ip), i=1,nflds)
            enddo
          elseif(nid.eq.ipass) then
            call csend(ipass,fieldout,len,0,nid)
          endif

        else  !ipass.eq.npass

          if(nid.eq.0) then
            do ip = 1,il
              write(50,'(1p20E14.6)'),
c     &         (pts(i,ip), i=1,ndim),
     &         (fieldout(i,ip), i=1,nflds)
            enddo
          endif

        endif
      enddo
      call nekgsync
      if(nid.eq.0)then
        close(unit=50)
        iFileNum=iFileNum+1
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine ps_hpts_out_fld(prefix,fieldout,nflds,
     $                             nfldm,npoints,nbuff)

c     output .fld file 

      include 'SIZE'
      include 'TOTAL'
      include 'RESTART'
C
C     Work arrays and temporary arrays
C
      common /scrcg/ pm1 (lx1,ly1,lz1,lelv)
      common/hpts_to_elm/NELGH,NXH,NYH,NZH
c
c     note, this usage of CTMP1 will be less than elsewhere if NELT ~> 3.
      parameter (lxyz=lx1*ly1*lz1)
      parameter (lpsc9=ldimt1+9)
      common /cbuff1/ tbuf(lxyz,lpsc9)
      real*4         tbuf

      real*4         test_pattern

      character*3    prefix
      real buf(nfldm,nbuff),fieldout(nfldm,nbuff)!from hpts_out
      character*1    fhdfle1(132)
      character*132   fhdfle
      equivalence   (fhdfle,fhdfle1)
      character*1    fldfile2(120)
      integer        fldfilei( 60)
      equivalence   (fldfilei,fldfile2)

      character*1 excode(30)
      character*10 frmat

      common /nopenf/ nopen(99)

      common /rdump/ ntdump
      data ndumps / 0 /

      logical ifxyo_s
      integer hxyz,ncount
      len = wdsize*nfldm*nbuff!from hpts_out
      hxyz=nxh*nyh*nzh

      npass = npoints/nbuff + 1
      il = mod(npoints,nbuff)
      if(il.eq.0) then
         il = nbuff
         npass = npass-1
      endif

c  Write to logfile that you're outputting data
      if(nio.eq.0) then 
        WRITE(6,1001) istep,time
 1001   FORMAT(/,i9,1pe12.4,' Write checkpoint:')
      endif
      call nekgsync()      

c  Check for file type
c  If filetype =6 then use multi-file-output
      p66 = abs(param(66))
c      if (p66.eq.6) then
c         call mfo_outfld(prefix)
c         call nekgsync                ! avoid race condition w/ outfld
c         return
c      endif

      ifxyo_s = ifxyo              ! Save ifxyo

c  Check given prefix against the database of prefixes
      iprefix = i_find_prefix(prefix,99)

      ierr = 0
      if (nid.eq.0) then

c       Open new file for each dump on /cfs
        nopen(iprefix)=nopen(iprefix)+1

        if (prefix.eq.'   '.and.nopen(iprefix).eq.1) ifxyo = .true. ! 1st file

        if (prefix.eq.'rst'.and.max_rst.gt.0) 
     $         nopen(iprefix) = mod1(nopen(iprefix),max_rst) ! restart

        call file2(nopen(iprefix),prefix)
c     if file type is 0 or negative then open using statement for ASCII
        if (p66.lt.1.0) then
           open(unit=24,file=fldfle,form='formatted',status='unknown')
        else
c     open binary file
           call  izero    (fldfilei,33)
           len = ltrunc   (fldfle,131)
           call chcopy    (fldfile2,fldfle,len)
           call byte_open (fldfile2,ierr)
c          write header as character string
           call blank(fhdfle,132)
        endif
      endif
c    broadcast if you are dumping the grid
      call bcast(ifxyo,lsize)
c    check to see if there was an error when byte_open was called
      if(p66.ge.1.0)
     $   call err_chk(ierr,'Error opening file in outfld. Abort. $')

C     Figure out what goes in EXCODE (header)
      CALL BLANK(EXCODE,30)
      NDUMPS=NDUMPS+1
      i=1
      if (mod(p66,1.0).eq.0.0) then !old header format
         IF(IFXYO) then
            EXCODE(1)='X'
            EXCODE(2)=' '
            EXCODE(3)='Y'
            EXCODE(4)=' '
            i = 5
            IF(IF3D) THEN
              EXCODE(i)  ='Z'
              EXCODE(i+1)=' '
              i = i + 2
            ENDIF
         ENDIF
         IF(IFVO) then
            EXCODE(i)  ='U'
            EXCODE(i+1)=' '
            i = i + 2
         ENDIF
         IF(IFPO) THEN
           EXCODE(i)='P'
           EXCODE(i+1)=' '
           i = i + 2
         ENDIF
         IF(IFTO) THEN
           EXCODE(i)='T '
           EXCODE(i+1)=' '
           i = i + 1
         ENDIF
         do iip=1,ldimt1
            if (ifpsco(iip)) then
              write(excode(iip+I)  ,'(i1)') iip
              write(excode(iip+I+1),'(a1)') ' '
              i = i + 1 
            endif
         enddo
      else
         !new header format
         IF (IFXYO) THEN !dumping grid information
            EXCODE(i)='X'
            i = i + 1
         ENDIF
         IF (IFVO) THEN !dumping velocity information
            EXCODE(i)='U'
            i = i + 1
         ENDIF
         IF (IFPO) THEN !dumping pressure information
            EXCODE(i)='P'
            i = i + 1
         ENDIF
         IF (IFTO) THEN !dumping Temperature information
            EXCODE(i)='T'
            i = i + 1
         ENDIF
         IF (LDIMT.GT.1) THEN !dumping passive scalar information
            NPSCALO = 0
            do k = 1,ldimt-1
              if(ifpsco(k)) NPSCALO = NPSCALO + 1
            enddo
            IF (NPSCALO.GT.0) THEN
               EXCODE(i) = 'S'
               WRITE(EXCODE(i+1),'(I1)') NPSCALO/10
               WRITE(EXCODE(i+2),'(I1)') NPSCALO-(NPSCALO/10)*10
            ENDIF
         ENDIF
      endif
     
c^^^^^^^^^^^^^^^^^^^^^^No changes necessary^^^^^^^^^^^^^^^^^^^^^^^^^^^^
c !!!!!!!!!!!!!!!!!!!!!Begin Changes!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c
c     Begining from hpts_out set up amount to pass
      npass = npoints/nbuff + 1
      il = mod(npoints,nbuff)
      if(il.eq.0) then
         il = nbuff
         npass = npass-1
      endif
C     Dump header based on phil files
      ierr = 0
      if (nid.eq.0) call ps_dump_header(excode,p66,ierr)
      call err_chk(ierr,'Error dumping header in outfld. Abort. $')

c   Get number of fields to write to file (xyzuvwpTt1 etc)
      call get_id(id)

      ierr = 0
      ncount=1
c     Dump out hpts in terms of elements
      do ipass = 1,npass

        if(ipass.lt.npass) then
          if(nid.eq.0) then
            call crecv(ipass,buf,len)
            do ip = 1,nbuff
              do i=1,nflds
                tbuf(ncount,i)=buf(i,ip) 
              enddo
              ncount=ncount+1
              if(ncount-1.eq.hxyz)then
                 call ps_out_buff(id,p66,ierr)
                 ncount=1
              endif
            enddo
          elseif(nid.eq.ipass) then
            call csend(ipass,fieldout,len,0,nid)
          endif

        else  !ipass.eq.npass

          if(nid.eq.0) then
            do ip = 1,il
              do i=1,nflds
               tbuf(ncount,i)=fieldout(i,ip)
              enddo
              ncount=ncount+1
              if(ncount-1.eq.hxyz)then
                 call ps_out_buff(id,p66,ierr)
                 ncount=1
              endif
            enddo
          endif

        endif
      enddo
      call nekgsync

      call err_chk(ierr,'Error writing file in outfld. Abort. $')

      ifxyo = ifxyo_s           ! restore ifxyo

      if (nid.eq.0) call close_fld(p66,ierr)
      call err_chk(ierr,'Error closing file in outfld. Abort. $')

      return
      end
c-----------------------------------------------------------------------
      subroutine ps_dump_header(excodein,p66,ierr)

      include 'SIZE'
      include 'TOTAL'
      common/hpts_to_elm/NELGH,NXH,NYH,NZH

      character*30  excodein

      character*30 excode
      character*1  excode1(30)
      equivalence (excode,excode1) 

      real*4         test_pattern

      character*1 fhdfle1(132)
      character*132 fhdfle
      equivalence (fhdfle,fhdfle1)

      write(excode,'(A30)') excodein

      ikstep = istep
      do ik=1,10
         if (ikstep.gt.9999) ikstep = ikstep/10
      enddo

      call blank(fhdfle,132)

c       write(6,111)               !       print on screen
c     $     nelgt,nx1,ny1,nz1,time,istep,excode
c
      if (mod(p66,1.0).eq.0.0) then !       old header format
         if (p66.lt.1.0) then       !ASCII
           if(nelgh.lt.10000) then
            WRITE(24,'(4i4,1pe14.7,I5,1X,30A1,1X,A12)')
     $           NELGH,NXH,NYH,NZH,TIME,ikstep,(EXCODE1(I),I=1,30),
     $           'NELT,NX,NY,N'
           else
            WRITE(24,'(i10,3i4,1pe18.9,I9,1X,30A1,1X,A12)')
     $           NELGH,NXH,NYH,NZH,TIME,ikstep,(EXCODE1(I),I=1,30),
     $           'NELT,NX,NY,N'
           endif
         else                       !Binary
            if (nelgh.lt.10000) then
               WRITE(fhdfle,'(4I4,1pe14.7,I5,1X,30A1,1X,A12)')
     $              NELGH,NXH,NYH,NZH,TIME,ikstep,(EXCODE1(I),I=1,30),
     $              ' 4 NELT,NX,NY,N'
            else
               write(fhdfle,'(i10,3i4,1P1e18.9,i9,1x,30a1)')
     $         nelgh,nxh,nyh,nzh,time,istep,(excode1(i),i=1,30)
            endif
            call byte_write(fhdfle,20,ierr)
         endif
      else                        !       new header format
         if (p66.eq.0.1) then
            write(24,111)
     $           nelgh,nxh,nyh,nzh,time,istep,excode
        else       
             write(fhdfle,111)
     $            nelgh,nxh,nyh,nzh,time,istep,excode
             call byte_write(fhdfle,20,ierr)
        endif
 111    FORMAT(i10,1x,i2,1x,i2,1x,i2,1x,1P1e18.9,1x,i9,1x,a)
      endif

      if(ierr.ne.0) return

      CDRROR=0.0
      if (p66.LT.1.0) then       !       formatted i/o
         WRITE(24,'(6G11.4)')(CDRROR,I=1,NELGH)   ! dummy 
      else
C       write byte-ordering test pattern to byte file...
        test_pattern = 6.54321
        call byte_write(test_pattern,1,ierr)
      endif

      return
      end
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
      subroutine ps_out_buff(id,p66,ierr)

      include 'SIZE'
      include 'TOTAL'

      parameter (lpsc9=ldimt1+9)
      parameter (lxyz=lx1*ly1*lz1)

      common /cbuff1/ tbuf(lxyz,lpsc9)
      common/hpts_to_elm/NELGH,NXH,NYH,NZH
      real*4         tbuf

      character*11 frmat
      nxyz=NXH*NYH*NZH
      
      call blank(frmat,11)
      if (id.le.9) then
         WRITE(FRMAT,1801) ID
 1801    FORMAT('(1p',I1,'e14.6)')
      else
         WRITE(FRMAT,1802) ID
 1802    FORMAT('(1p',I2,'e14.6)')
      endif

      if (p66.lt.1.0) then
C       formatted i/o
        WRITE(24,FRMAT)
     $      ((TBUF(I,II),II=1,ID),I=1,nxyz)
      else
C        C binary i/o
         do ii=1,id
c            call byte_reverse(tbuf,id,ierr)
            call byte_write(tbuf(1,ii),nxyz,ierr)
            if(ierr.ne.0) goto 101
         enddo
      endif
 101  continue

      return
      end
c-----------------------------------------------------------------------
