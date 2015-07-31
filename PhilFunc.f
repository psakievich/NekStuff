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
      subroutine ps_hpts
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

      parameter(nfldm=ldim+ldimt+1)


      common /c_hptsr/ pts      (ldim,lhis)
     $               , fieldout (nfldm,lhis)
     $               , dist     (lhis)
     $               , rst      (lhis*ldim)


      common /c_hptsi/ rcode(lhis),elid(lhis),proc(lhis)

      common /scrcg/  pm1 (lx1,ly1,lz1,lelv) ! mapped pressure
      common /outtmp/ wrk (lx1*ly1*lz1*lelt,nfldm)

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
        call ps_hpts_in(pts,npts,npoints)
        call intpts_setup(-1.0,inth_hpts) ! use default tolerance
      endif


      call prepost_map(0)  ! maps axisymm and pressure

      ! pack working array
      nflds = 0
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

      ! evaluate input field at given points
      do ifld = 1,nflds
         call findpts_eval(inth_hpts,fieldout(ifld,1),nfldm,
     &                     rcode,1,
     &                     proc,1,
     &                     elid,1,
     &                     rst,ndim,npts,
     &                     wrk(1,ifld))
      enddo
      ! write interpolation results to file
      call ps_hpts_out(pts,fieldout,nflds,nfldm,npoints,nbuff)

      call prepost_map(1)  ! maps back axisymm arrays

      if(nio.eq.0) write(6,*) 'done :: dump history points'

      return
      end
c-----------------------------------------------------------------------
      subroutine ps_buffer_in(buffer,npp,npoints,nbuf)
        
      include 'SIZE'
      include 'PARALLEL'

      real    buffer(ldim,nbuf)  

      ierr = 0
      if(nid.eq.0) then
        write(6,*) 'reading hpts.in'
        open(50,file='hpts.in',status='old',err=100)
        read(50,*,err=100) npoints
        goto 101
 100    ierr = 1
 101    continue
      endif
      ierr=iglsum(ierr,1)
      if(ierr.gt.0) then
        write(6,*) 'Cannot open hpts.in in subroutine hpts()'
        call exitt
      endif
      
      call bcast(npoints,isize)
      if(npoints.gt.lhis*np) then
        if(nid.eq.0) write(6,*) 'ABORT: Too many pts to read in hpts()!'
        call exitt
      endif
      if(nid.eq.0) write(6,*) 'found ', npoints, ' points'


      npass =  npoints/nbuf +1  !number of passes to cover all pts
      n0    =  mod(npoints,nbuf)!remainder 
      if(n0.eq.0) then
         npass = npass-1
         n0    = nbuf
      endif

      len = wdsize*ndim*nbuf
      if (nid.gt.0.and.nid.lt.npass) msg_id=irecv(nid,buffer,len)
      call nekgsync
      
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
      real    pts(ldim,npts)

      if (lt2.gt.npts) then

         call ps_buffer_in(xyz,npp,npoints,lt2)
         if(npoints.gt.np*npts) then
           if(nid.eq.0)write(6,*)'ABORT in hpts(): npoints > NP*lhis!!' 
           if(nid.eq.0)write(6,*)'Change SIZE: ',np,npts,npoints
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
      subroutine ps_hpts_out(pts,fieldout,nflds,nfldm,npoints,nbuff)
c    ********************************************
c    *** MODIFIED VERSION OF HPTS_OUT IN REPO****
c    ********************************************

      include 'SIZE'
      include 'TOTAL'
      real pts(ldim,lhis)
      real buf(nfldm,nbuff),fieldout(nfldm,nbuff)
      character*80 filename
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

      do ipass = 1,npass
        if(nid.eq.0)then
          write(filename,"('cyl_coor',I0,'.dat')")iFileNum
          open(unit=50,file=filename)!,status='new')
c          write(50,'(A)') 
c     &       '# time  vx  vy  [vz]  pr  T  PS1  PS2  ...'
        endif
        call nekgsync

        if(ipass.lt.npass) then
          if(nid.eq.0) then
            call crecv(ipass,buf,len)
            do ip = 1,nbuff
              write(50,'(1p20E15.7)'),
     &         (pts(i,ip), i=1,ndim),
     &         (buf(i,ip), i=1,nflds)
            enddo
          elseif(nid.eq.ipass) then
            call csend(ipass,fieldout,len,0,nid)
          endif

        else  !ipass.eq.npass

          if(nid.eq.0) then
            do ip = 1,il
              write(50,'(1p20E15.7)'),
     &         (pts(i,ip), i=1,ndim),
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
