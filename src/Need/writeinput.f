      subroutine writeinput(unitno,lstring,rstring,tag)
c--- This routine echoes the line in the input file specified by the
c--- input parameter 'tag'; if 'tag' is set to 'WRITEALL' then all
c--- of the input file is echoed.
c--- Output is written to the unit 'unitno', with lines bracketed
c--- by the strings 'lstring' and 'rstring'
      implicit none
      include 'constants.f'
      include 'order.f'
      include 'maxwt.f'
      include 'masses.f'
      include 'zerowidth.f'
      include 'flags.f'
      include 'clustering.f'
      include 'anomcoup.f'
      include 'gridinfo.f'
      include 'limits.f'
      include 'jetcuts.f'
      include 'leptcuts.f'
      include 'lhapdf.f'
      include 'pdlabel.f'
      include 'removebr.f'
      include 'dynamicscale.f'
      include 'verbose.f'
      include 'debug.f'
      include 'virtonly.f'
      include 'realonly.f'
      include 'noglue.f'
      include 'realwt.f'
      include 'cutoff.f'
      include 'initialscales.f'
      include 'outputoptions.f'
      include 'outputflags.f'
      include 'part.f'
      include 'runstring.f'
      character*(*) tag,lstring,rstring
      character*72 f94,f95,f96,f97,f98,f99
      logical dryrun,makecuts,writeall,spira
      integer unitno, nmin,nmax
      integer nproc,ih1,ih2,itmx1,itmx2,ncall1,ncall2,origij
      integer NPTYPE,NGROUP,NSET
      double precision sqrts,rtsmin,Rcut
 
      common/spira/spira
      common/nmin/nmin
      common/nmax/nmax
      common/rtsmin/rtsmin


      common/nproc/nproc
      common/energy/sqrts
      common/density/ih1,ih2
      common/iterat/itmx1,ncall1,itmx2,ncall2
      common/dryrun/dryrun
      
      common/pdflib/NPTYPE,NGROUP,NSET
      
      common/Rcut/Rcut
      common/makecuts/makecuts

      common/origij/origij
      
c--- f94 integer.XXYY format      
      f94='('''//lstring//''',11x,i4,''.'',2a2,12x,''['',a,'']'','''
     & //rstring//''')' 
c--- f95 2xfloating point format
      f95='('''//lstring//''',f8.3,'','',f8.3,16x,''['',a,'']'','''
     & //rstring//''')' 
c--- f96 character format            
      f96='('''//lstring//''',a20,12x,''['',a,'']'','''//rstring//''')' 
c--- f97 integer format      
      f97='('''//lstring//''',i20,12x,''['',a,'']'','''//rstring//''')' 
c--- f98 logical format      
      f98='('''//lstring//''',L20,12x,''['',a,'']'','''//rstring//''')' 
c--- f99 floating point format
      f99='('''//lstring//''',f20.4,12x,''['',a,'']'','''
     & //rstring//''')' 
    
      writeall=.false.
      if (tag .eq. 'WRITEALL') writeall=.true.
      
      if (writeall) then
      write(unitno,*) lstring//' Run corresponds to this input file)'
      write(unitno,*)
      write(unitno,*)
     . lstring//' [Flags to specify the mode in which NNLOCAL is run] )'
      endif
      
      if (writeall) then
      write(unitno,*)
      write(unitno,*) lstring//
     & ' [General options to specify the process and execution] )'
      endif
      if ((tag .eq. 'nproc') .or. (writeall)) then
         write(unitno,fmt=f97) nproc,'nproc'
      endif
      if ((tag .eq. 'order') .or. (writeall)) then
      write(unitno,fmt=f97) order,'order'
      endif
      if ((tag .eq. 'part') .or. (writeall)) then
      write(unitno,fmt=f96) part,'part'
      endif
      if ((tag .eq. 'runstring') .or. (writeall)) then
      write(unitno,fmt=f96) runstring,'runstring'
      endif
      if ((tag .eq. 'sqrts') .or. (writeall)) then
      write(unitno,fmt=f99) sqrts,'sqrts'
      endif
      if ((tag .eq. 'ih1') .or. (writeall)) then
      write(unitno,fmt=f97) ih1,'ih1'
      endif
      if ((tag .eq. 'ih2') .or. (writeall)) then
      write(unitno,fmt=f97) ih2,'ih2'
      endif
      if ((tag .eq. 'hmass') .or. (writeall)) then
      write(unitno,fmt=f99) hmass,'hmass'
      endif
      if ((tag .eq. 'scale') .or. (writeall)) then
         write(unitno,fmt=f99) initscale,'scale'
      endif
      if ((tag .eq. 'facscale') .or. (writeall)) then
         write(unitno,fmt=f99) initfacscale,'facscale'
      endif
      if ((tag .eq. 'dynamicscale') .or. writeall) then
      write(unitno,fmt=f96) dynstring,'dynamicscale'
      endif
      if ((tag .eq. 'zerowidth') .or. (writeall)) then
      write(unitno,fmt=f98) zerowidth,'zerowidth'
      endif
      if ((tag .eq. 'removebr') .or. (writeall)) then
      write(unitno,fmt=f98) removebr,'removebr'
      endif
      if ((tag .eq. 'itmx1') .or. (writeall)) then
      write(unitno,fmt=f97) itmx1,'itmx1'
      endif
      if ((tag .eq. 'ncall1') .or. (writeall)) then
      write(unitno,fmt=f97) ncall1,'ncall1'
      endif
      if ((tag .eq. 'itmx2') .or. (writeall)) then
      write(unitno,fmt=f97) itmx2,'itmx2'
      endif
      if ((tag .eq. 'ncall2') .or. (writeall)) then
      write(unitno,fmt=f97) ncall2,'ncall2'
      endif
      if ((tag .eq. 'ij') .or. (writeall)) then
      write(unitno,fmt=f97) origij,'ij'
      endif
      if ((tag .eq. 'dryrun') .or. (writeall)) then
      write(unitno,fmt=f98) dryrun,'dryrun'
      endif
      if (writeall) then
      write(unitno,*)
      write(unitno,*) 
     . lstring//' [Heavy quark masses] )'
      endif
      if ((tag .eq. 'top mass') .or. (writeall)) then
      write(unitno,fmt=f99) mt,'top mass'
      endif
      if ((tag .eq. 'bottom mass') .or. (writeall)) then
      write(unitno,fmt=f99) mb,'bottom mass'
      endif
      if ((tag .eq. 'charm mass') .or. (writeall)) then
      write(unitno,fmt=f99) mc,'charm mass'
      endif

      if (writeall) then
      write(unitno,*)
      write(unitno,*) 
     . lstring//' [Pdf selection] )'
      endif
      if ((tag .eq. 'LHAPDF group') .or. (writeall)) then
      write(unitno,fmt=f96) PDFname,'LHAPDF group'
      endif
      if ((tag .eq. 'LHAPDF set') .or. (writeall)) then
      write(unitno,fmt=f97) PDFmember,'LHAPDF set'
      endif

      if (writeall) then
      write(unitno,*)
      write(unitno,*)
     . lstring//' [Jet definition and event cuts] )'
      endif
      if ((tag .eq. 'm34min') .or. (writeall)) then
      write(unitno,fmt=f99) dsqrt(wsqmin),'m34min'
      endif
      if ((tag .eq. 'm34max') .or. (writeall)) then
      write(unitno,fmt=f99) dsqrt(wsqmax),'m34max'
      endif
      if ((tag .eq. 'm56min') .or. (writeall)) then
      write(unitno,fmt=f99) dsqrt(bbsqmin),'m56min'
      endif
      if ((tag .eq. 'm56max') .or. (writeall)) then
      write(unitno,fmt=f99) dsqrt(bbsqmax),'m56max'
      endif
      if ((tag .eq. 'inclusive') .or. (writeall)) then
      write(unitno,fmt=f98) inclusive,'inclusive'
      endif
      if ((tag .eq. 'algorithm') .or. (writeall)) then
      write(unitno,fmt=f96) algorithm,'algorithm'
      endif
      if ((tag .eq. 'ptjetmin') .or. (writeall)) then
      write(unitno,fmt=f99) ptjetmin,'ptjetmin'
      endif
      if ((tag .eq. 'etajetmin') .or. (writeall)) then
      write(unitno,fmt=f99) etajetmin,'etajetmin'
      endif
      if ((tag .eq. 'etajetmax') .or. (writeall)) then
      write(unitno,fmt=f99) etajetmax,'etajetmax'
      endif
      if ((tag .eq. 'Rcut') .or. (writeall)) then
      write(unitno,fmt=f99) Rcut,'Rcut'
      endif
      if ((tag .eq. 'makecuts') .or. (writeall)) then
      write(unitno,fmt=f98) makecuts,'makecuts'
      endif
      if ((tag .eq. 'leptpt') .or. (writeall)) then
      write(unitno,fmt=f99) leptpt,'leptpt'
      endif
      if ((tag .eq. 'leptrap') .or. (writeall)) then
      write(unitno,fmt=f99) leptrap,'leptrap'
      endif
      if ((tag .eq. 'leptveto') .or. (writeall)) then
      write(unitno,fmt=f95) leptveto1min,leptveto1max,'leptveto'
      endif
      if ((tag .eq. 'misspt') .or. (writeall)) then
      write(unitno,fmt=f99) misspt,'misspt'
      endif
      if ((tag .eq. 'leptpt2') .or. (writeall)) then
      write(unitno,fmt=f99) leptpt2,'leptpt2'
      endif
      if ((tag .eq. 'leptrap2') .or. (writeall)) then
      write(unitno,fmt=f99) leptrap2,'leptrap2'
      endif
      if ((tag .eq. 'leptveto2') .or. (writeall)) then
      write(unitno,fmt=f95) leptveto2min,leptveto2max,'leptveto2'
      endif
      if ((tag .eq. 'mtrans34cut') .or. (writeall)) then
      write(unitno,fmt=f99) mtrans34cut,'mtrans34cut'
      endif
      if ((tag .eq. 'Rjlmin') .or. (writeall)) then
      write(unitno,fmt=f99) Rjlmin,'Rjlmin'
      endif
      if ((tag .eq. 'Rllmin') .or. (writeall)) then
      write(unitno,fmt=f99) Rllmin,'Rllmin'
      endif
      if ((tag .eq. 'delyjjmin') .or. (writeall)) then
      write(unitno,fmt=f99) delyjjmin,'delyjjmin'
      endif
      if ((tag .eq. 'jetsopphem') .or. (writeall)) then
      write(unitno,fmt=f98) jetsopphem,'jetsopphem'
      endif
      if ((tag .eq. 'lbjscheme') .or. (writeall)) then
      write(unitno,fmt=f97) lbjscheme,'lbjscheme'
      endif
      if ((tag .eq. 'ptbjetmin') .or. (writeall)) then
      write(unitno,fmt=f99) ptbjetmin,'ptbjetmin'
      endif
      if ((tag .eq. 'etabjetmax') .or. (writeall)) then
      write(unitno,fmt=f99) etabjetmax,'etabjetmax'
      endif

      if (writeall) then
      write(unitno,*)
      write(unitno,*) 
     . lstring//' [How to resume/save a run] )'
      endif
      if ((tag .eq. 'readin') .or. (writeall)) then
      write(unitno,fmt=f98) readin,'readin'
      endif
      if ((tag .eq. 'writeout') .or. (writeall)) then
      write(unitno,fmt=f98) writeout,'writeout'
      endif

      if (writeall) then
      write(unitno,*)
      write(unitno,*) lstring//
     & ' [Technical parameters that should not normally be changed]'
      endif
      if ((tag .eq. 'debug') .or. (writeall)) then
      write(unitno,fmt=f98) debug,'debug'
      endif
      if ((tag .eq. 'verbose') .or. (writeall)) then
      write(unitno,fmt=f98) verbose,'verbose'
      endif
      if ((tag .eq. 'virtonly') .or. (writeall)) then
      write(unitno,fmt=f98) virtonly,'virtonly'
      endif
      if ((tag .eq. 'realonly') .or. (writeall)) then
      write(unitno,fmt=f98) realonly,'realonly'
      endif
      if ((tag .eq. 'spira') .or. (writeall)) then
      write(unitno,fmt=f98) spira,'spira'
      endif
      if ((tag .eq. 'noglue') .or. (writeall)) then
      write(unitno,fmt=f98) noglue,'noglue'
      endif
      if ((tag .eq. 'ggonly') .or. (writeall)) then
      write(unitno,fmt=f98) ggonly,'ggonly'
      endif
      if ((tag .eq. 'gqonly') .or. (writeall)) then
      write(unitno,fmt=f98) gqonly,'gqonly'
      endif
      if ((tag .eq. 'omitgg') .or. (writeall)) then
      write(unitno,fmt=f98) omitgg,'omitgg'
      endif
      if ((tag .eq. 'nmin') .or. (writeall)) then
      write(unitno,fmt=f97) nmin,'nmin'
      endif
      if ((tag .eq. 'nmax') .or. (writeall)) then
      write(unitno,fmt=f97) nmax,'nmax'
      endif
      if ((tag .eq. 'clustering') .or. (writeall)) then
      write(unitno,fmt=f98) clustering,'clustering'
      endif
      if ((tag .eq. 'realwt') .or. (writeall)) then
      write(unitno,fmt=f98) realwt,'realwt'
      endif
      if ((tag .eq. 'rtsmin') .or. (writeall)) then
      write(unitno,fmt=f99) rtsmin,'rtsmin'
      endif
      if ((tag .eq. 'cutoff') .or. (writeall)) then
      write(unitno,fmt=f99) cutoff,'cutoff'
      endif
      
      if (writeall) then
      write(unitno,*)
      endif
      
      return

c--- 96 character format      
c   96 format(' (',a20,12x,'[',a,']',' )')  
c--- 97 integer format      
c   97 format(' (',i20,12x,'[',a,']',' )')  
c--- 98 logical format      
c   98 format(' (',L20,12x,'[',a,']',' )')  
c--- 99 floating point format
c   99 format(' (',f20.4,12x,'[',a,']',' )')  
      
      end
      
