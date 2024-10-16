      subroutine hp_iopmat2(p,x1,x2)
      implicit none
      include 'hp_types.h'
      include 'hp_constants.h'
      include 'hp_Ioperators2.f'
      include 'hp_APinc.h'
      include 'order.f'
      include 'scale.f'
      include 'facscale.f'
      include 'agq.f'
      real(ki) b0
      integer i,j,k
      real(ki) xl12,xlfs,xlrf,p(mxpart,4)
      double precision x1,x2,dot
      double precision xl12_lowprec,xlfs_lowprec,xlrf_lowprec
      real(ki) xa,xb
      real(ki) resab(2,2,-2:2),resba(2,2,-2:2)
      real(ki) tmpab(2,2,-4:0),tmpba(2,2,-4:0)

      real(ki) res2(-4:0),res3(-4:0),res(2,2,-4:0)
      real(ki) tmp1(2,2,-2:2),tmp2(2,2,-4:0),tmp(2,2,-4:0),tmp3(0:2,0:2,-4:0)
      real(ki) tmpab1(2,2,-2:2),tmpab2(2,2,-4:0),tmpab3(0:2,0:2,-4:0)
      real(ki) tmpba1(2,2,-2:2),tmpba2(2,2,-4:0),tmpba3(0:2,0:2,-4:0)

      real(ki) ABa2carsggg(2,2,-4:0),BAa2carsggg(2,2,-4:0)
      real(ki) ABa2carbsgggg(2,2,-4:0),BAa2carbsgggg(2,2,-4:0)
      real(ki) ABa2srsgg(0:2,0:2,-4:0),BAa2srsgg(0:2,0:2,-4:0)
      real(ki) ABa2carssrsgg(0:2,0:2,-4:0),BAa2carssrsgg(0:2,0:2,-4:0)
      real(ki) ABa2carbssrsgg(0:2,0:2,-4:0),BAa2carbssrsgg(0:2,0:2,-4:0)

      real(ki) ABa12carsssggg(0:2,0:2,-4:0),BAa12carsssggg(0:2,0:2,-4:0)   
      real(ki) ABa12carscrsssggg(0:2,0:2,-4:0),BAa12carscrsssggg(0:2,0:2,-4:0)
      real(ki) ABa12carscasssggg(0:2,0:2,-4:0),BAa12carscasssggg(0:2,0:2,-4:0)
      real(ki) ABa12carscrsggg(2,2,-4:0),BAa12carscrsggg(2,2,-4:0)   
      real(ki) ABa12carscasggg(2,2,-4:0),BAa12carscasggg(2,2,-4:0)   
      real(ki) ABa12casbrcasgggg(2,2,-4:0),BAa12casbrcasgggg(2,2,-4:0)
      real(ki) ABa12srscrsgg(0:2,0:2,-4:0),BAa12srscrsgg(0:2,0:2,-4:0)   
      real(ki) ABa12carssrscrsgg(0:2,0:2,-4:0),BAa12carssrscrsgg(0:2,0:2,-4:0)  
      real(ki) ABa12carssrsssgg(0:2,0:2,-4:0),BAa12carssrsssgg(0:2,0:2,-4:0)    
      real(ki) ABa12srscasssgg(0:2,0:2,-4:0),BAa12srscasssgg(0:2,0:2,-4:0)      
      real(ki) ABa12srsssgg(0:2,0:2,-4:0),BAa12srsssgg(0:2,0:2,-4:0)
      real(ki) ABa12carssrscasssgg(0:2,0:2,-4:0),BAa12carssrscasssgg(0:2,0:2,-4:0)
      
      real(ki) ABa1a1crssscarggg(0:2,0:2,-4:0),BAa1a1crssscarggg(0:2,0:2,-4:0) 
      real(ki) ABa1a1sscarggg(0:2,0:2,-4:0),BAa1a1sscarggg(0:2,0:2,-4:0)    
      real(ki) ABa1a1carcasggg(2,2,-4:0),BAa1a1carcasggg(2,2,-4:0)   
      real(ki) ABa1a1cbrcasgggg(2,2,-4:0),BAa1a1cbrcasgggg(2,2,-4:0)  
      real(ki) ABa1a1crscasggg(2,2,-4:0),BAa1a1crscasggg(2,2,-4:0)   
      real(ki) ABa1a1cassscarsragg(0:2,0:2,-4:0),BAa1a1cassscarsragg(0:2,0:2,-4:0)
      real(ki) ABa1a1cassssragg(0:2,0:2,-4:0),BAa1a1cassssragg(0:2,0:2,-4:0)  
      real(ki) ABa1a1sscarsragg(0:2,0:2,-4:0),BAa1a1sscarsragg(0:2,0:2,-4:0)        
      real(ki) ABa1a1sssrgg(0:2,0:2,-4:0),BAa1a1sssrgg(0:2,0:2,-4:0)            
      
      real(ki) ABrva1ncar1gg(2,2,-4:0),BArva1ncar1gg(2,2,-4:0)
      real(ki) ABrva1sr1g(0:2,0:2,-4:0),BArva1sr1g(0:2,0:2,-4:0)
      real(ki) ABrva1carsr1g(0:2,0:2,-4:0),BArva1carsr1g(0:2,0:2,-4:0)

      real(ki) ABa1ncar0gg(2,2,-2:2),BAa1ncar0gg(2,2,-2:2)
c      real(ki) ABa1sr0g(2,2,-2:2),BAa1sr0g(2,2,-2:2)
c      real(ki) ABa1carsr0g(2,2,-2:2),BAa1carsr0g(2,2,-2:2)

      real(ki) ABa1c1cargggg(2,2,-2:2),BAa1c1cargggg(2,2,-2:2)
      real(ki) ABa1c1cbrgggg(2,2,-2:2),BAa1c1cbrgggg(2,2,-2:2)

      real(ki) ABPgg0(2,2), ABPgq0(2,2), ABPqg0(2,2), ABPqq0(2,2)
      real(ki) BAPgg0(2,2), BAPgq0(2,2), BAPqg0(2,2), BAPqq0(2,2)
      real(ki) ABPgg1(2,2), ABPgq1(2,2), ABPqg1(2,2), ABPqq1(2,2)
      real(ki) BAPgg1(2,2), BAPgq1(2,2), BAPqg1(2,2), BAPqq1(2,2)
      real(ki) ABP0ggxP0gg(2,2), ABP0ggxP0gq(2,2), ABP0ggxP0qg(2,2), ABP0ggxP0qq(2,2),
     1         ABP0gqxP0gg(2,2), ABP0gqxP0gq(2,2), ABP0gqxP0qg(2,2), ABP0gqxP0qq(2,2),
     2         ABP0qgxP0gg(2,2), ABP0qgxP0gq(2,2), ABP0qgxP0qg(2,2), ABP0qgxP0qq(2,2),
     3         ABP0qqxP0gg(2,2), ABP0qqxP0gq(2,2), ABP0qqxP0qg(2,2), ABP0qqxP0qq(2,2)
      real(ki) BAP0ggxP0gg(2,2), BAP0ggxP0gq(2,2), BAP0ggxP0qg(2,2), BAP0ggxP0qq(2,2),
     1         BAP0gqxP0gg(2,2), BAP0gqxP0gq(2,2), BAP0gqxP0qg(2,2), BAP0gqxP0qq(2,2),
     2         BAP0qgxP0gg(2,2), BAP0qgxP0gq(2,2), BAP0qgxP0qg(2,2), BAP0qgxP0qq(2,2),
     3         BAP0qqxP0gg(2,2), BAP0qqxP0gq(2,2), BAP0qqxP0qg(2,2), BAP0qqxP0qq(2,2)
      real(ki) ABP0ggP0gg(2,2), ABP0ggP0gq(2,2), ABP0ggP0qg(2,2), ABP0ggP0qq(2,2),
     1         ABP0gqP0gg(2,2), ABP0gqP0gq(2,2), ABP0gqP0qg(2,2), ABP0gqP0qq(2,2),
     2         ABP0qgP0gg(2,2), ABP0qgP0gq(2,2), ABP0qgP0qg(2,2), ABP0qgP0qq(2,2),
     3         ABP0qqP0gg(2,2), ABP0qqP0gq(2,2), ABP0qqP0qg(2,2), ABP0qqP0qq(2,2)
      real(ki) BAP0ggP0gg(2,2), BAP0ggP0gq(2,2), BAP0ggP0qg(2,2), BAP0ggP0qq(2,2),
     1         BAP0gqP0gg(2,2), BAP0gqP0gq(2,2), BAP0gqP0qg(2,2), BAP0gqP0qq(2,2),
     2         BAP0qgP0gg(2,2), BAP0qgP0gq(2,2), BAP0qgP0qg(2,2), BAP0qgP0qq(2,2),
     3         BAP0qqP0gg(2,2), BAP0qqP0gq(2,2), BAP0qqP0qg(2,2), BAP0qqP0qq(2,2)
 
      real(ki) tempfcn(2,2)


      character *60 FMT3
      integer eps
      
      FMT3="(A20,4(F30.17,','),F30.16,'};')"

      xa=real(x1,kind=ki)
      xb=real(x2,kind=ki)
c      xa=0.1_ki
c      xb=0.2_ki
      
c      write(6,*) 'x1 xa =',x1,xa
c      write(6,*) 'x2 xb =',x2,xb
c      pause
      
c---  xl12=Log[Q^2/muR^2]
      xl12_lowprec=log(two*dot(p,1,2)/musq)
      xl12=real(xl12_lowprec,kind=ki)
      
c---  xlfs=Log[Q^2/muF^2]
      xlfs_lowprec=log(two*dot(p,1,2)/facscale**2)
      xlfs=real(xlfs_lowprec,kind=ki)

c---  xlrf=Log[muR^2/muF^2]
      xlrf_lowprec=log(musq/facscale**2)
      xlrf=real(xlrf_lowprec,kind=ki)

c      ason2pi2=ason2pi**2
            
c---  Log[Q^2/muR^2]
c      xl12=zip

c---  xlfs=Log[Q^2/muF^2]
c      xlfs=zip

c---  xlrf=Log[muR^2/muF^2]
c      xlrf=zip
           
      b0=(xn*11.0_ki)/6.0_ki
           
c 10    continue

      ! (a,b,c,d) stand for incoming, hard, final1, final2
      do i=0,2
         do j=0,2
            do eps=-4,0
               I20op1eps(g,g,g,g,i,j,eps)=0.0_ki
               I20op2eps(g,g,g,g,i,j,eps)=0.0_ki
            enddo
         enddo
      enddo
      
      ! (a,b,c) stand for incoming, hard, final1, final2
      do i=1,2
         do j=1,2
            do eps=-2,2
               I10op1eps(g,g,g,i,j,eps)=0.0_ki
               I10op2eps(g,g,g,i,j,eps)=0.0_ki
            enddo
         enddo
      enddo


      
c---  A2 Integrated counterterms

      call hp_a2carsggg(xa,xb,ABa2carsggg)
      call hp_a2carbsgggg(xa,xb,ABa2carbsgggg)
      call hp_a2carsggg_bulk(xa,xb,res)
      ABa2carsggg(2,2,-4:0)=res(2,2,-4:0)
      call hp_a2carbsgggg_bulk(xa,xb,res)
      ABa2carbsgggg(2,2,-4:0)=res(2,2,-4:0)

      call hp_a2carsggg(xb,xa,BAa2carsggg)
      call hp_a2carbsgggg(xb,xa,BAa2carbsgggg)
      call hp_a2carsggg_bulk(xb,xa,res)
      BAa2carsggg(2,2,-4:0)=res(2,2,-4:0)
      call hp_a2carbsgggg_bulk(xb,xa,res)
      BAa2carbsgggg(2,2,-4:0)=res(2,2,-4:0)

      tmpab2=ABa2carsggg+ABa2carbsgggg
      tmpba2=BAa2carsggg+BAa2carbsgggg

      tmpab2=half*tmpab2
      tmpba2=half*tmpba2
            
      k=1
      call hp_addscale2(k,tmp2,tmpab2,xl12)
      call hp_fillIop2_1(g,g,g,g,k,tmp2)
      
      call hp_addscale2(k,tmp2,tmpba2,xl12)
      call hp_fillIop2_2(g,g,g,g,k,tmp2)

c---
      
      call hp_a2srsgg(CA,xa,ABa2srsgg)
      call hp_a2srsgg(CA,xb,BAa2srsgg)
      call hp_a2carssrsgg(CA,xa,ABa2carssrsgg)
      call hp_a2carssrsgg(CA,xb,BAa2carssrsgg)
      call hp_a2carbssrsgg(CA,CA,xa,ABa2carbssrsgg)
      call hp_a2carbssrsgg(CA,CA,xb,BAa2carbssrsgg)

      tmpab3=half*ABa2srsgg-ABa2carssrsgg-ABa2carbssrsgg
      tmpba3=half*BAa2srsgg-BAa2carssrsgg-BAa2carbssrsgg

      tmpab3=half*tmpab3
      tmpba3=half*tmpba3

      k=0
      call hp_addscale2(k,tmp3,tmpab3,xl12)
      call hp_fillIop2_1(g,g,g,g,k,tmp3)
      
      call hp_addscale2(k,tmp3,tmpba3,xl12)
      call hp_fillIop2_2(g,g,g,g,k,tmp3)

      
c$$$      write(6,*) 'xa=',xa
c$$$      write(6,*) 'xb=',xb
c$$$      write(6,*)
c$$$c      write(6,*) 'I10'
c$$$      write(6,*)
c$$$      write(6,FMT3) 'I10opeps[1,1]={',I10op1eps(g,g,g,1,1,:)+I10op2eps(g,g,g,1,1,:)
c$$$      write(6,*)
c$$$      write(6,FMT3) 'I10opeps[1,2]={',I10op1eps(g,g,g,1,2,:)+I10op2eps(g,g,g,2,1,:)
c$$$      write(6,*)
c$$$      write(6,FMT3) 'I10opeps[2,1]={',I10op1eps(g,g,g,2,1,:)+I10op2eps(g,g,g,1,2,:)
c$$$      write(6,*)
c$$$      write(6,FMT3) 'I10opeps[2,2]={',I10op1eps(g,g,g,2,2,:)+I10op2eps(g,g,g,2,2,:)
c$$$      write(6,*)
c$$$      write(6,*)
c$$$      write(6,*) 'I20'
c$$$      write(6,*)
c$$$      write(6,FMT3) 'I20opeps[1,0]={',I20op1eps(g,g,g,g,1,0,:)+I20op2eps(g,g,g,g,0,1,:)
c$$$      write(6,*)                   
c$$$      write(6,FMT3) 'I20opeps[0,1]={',I20op1eps(g,g,g,g,0,1,:)+I20op2eps(g,g,g,g,1,0,:)
c$$$      write(6,*)                   
c$$$      write(6,FMT3) 'I20opeps[0,0]={',I20op1eps(g,g,g,g,0,0,:)+I20op2eps(g,g,g,g,0,0,:)
c$$$c      write(6,*) 'I20op1eps[0,0]=',I20op1eps(g,g,g,g,0,0,:)
c$$$c      write(6,*) 'I20op2eps[0,0]=',I20op2eps(g,g,g,g,0,0,:)
c$$$      write(6,*)                   
c$$$      write(6,FMT3) 'I20opeps[1,1]={',I20op1eps(g,g,g,g,1,1,:)+I20op2eps(g,g,g,g,1,1,:)
c$$$      write(6,*)                   
c$$$      write(6,FMT3) 'I20opeps[1,2]={',I20op1eps(g,g,g,g,1,2,:)+I20op2eps(g,g,g,g,2,1,:)
c$$$      write(6,*)
c$$$      write(6,FMT3) 'I20opeps[2,1]={',I20op1eps(g,g,g,g,2,1,:)+I20op2eps(g,g,g,g,1,2,:)
c$$$      write(6,*)
c$$$      write(6,FMT3) 'I20opeps[2,2]={',I20op1eps(g,g,g,g,2,2,:)+I20op2eps(g,g,g,g,2,2,:)
c$$$      write(6,*) 'in Iopmat'
c$$$      pause



c      goto 20

      
      
c---  A12 Integrated counterterms

      call hp_a12carscrsggg(xa,xb,ABa12carscrsggg)
      call hp_a12carscrsggg(xb,xa,BAa12carscrsggg)
      call hp_a12carscasggg(xa,xb,ABa12carscasggg)
      call hp_a12carscasggg(xb,xa,BAa12carscasggg)
      call hp_a12casbrcasgggg(xa,xb,ABa12casbrcasgggg)
      call hp_a12casbrcasgggg(xb,xa,BAa12casbrcasgggg)

      tmpab2=
     1     two*ABa12carscasggg+two*ABa12casbrcasgggg+ABa12carscrsggg

      tmpba2=
     1     two*BAa12carscasggg+two*BAa12casbrcasgggg+BAa12carscrsggg

      tmpab2=-half*tmpab2
      tmpba2=-half*tmpba2

      k=1
      call hp_addscale2(k,tmp2,tmpab2,xl12)
      call hp_fillIop2_1(g,g,g,g,k,tmp2)

      call hp_addscale2(k,tmp2,tmpba2,xl12)
      call hp_fillIop2_2(g,g,g,g,k,tmp2)
      
c---

      call hp_a12srscrsgg(CA,xa,ABa12srscrsgg)
      call hp_a12srscrsgg(CA,xb,BAa12srscrsgg)
      call hp_a12carssrscrsgg(CA,xa,ABa12carssrscrsgg)
      call hp_a12carssrscrsgg(CA,xb,BAa12carssrscrsgg)
      call hp_a12carssrsssgg(CA,xa,ABa12carssrsssgg)
      call hp_a12carssrsssgg(CA,xb,BAa12carssrsssgg)
      call hp_a12srsssgg(CA,xa,ABa12srsssgg)
      call hp_a12srsssgg(CA,xb,BAa12srsssgg)
      call hp_a12carssrscasssgg(CA,xa,ABa12carssrscasssgg)
      call hp_a12carssrscasssgg(CA,xb,BAa12carssrscasssgg)
      call hp_a12srscasssgg(CA,xa,ABa12srscasssgg)
      call hp_a12srscasssgg(CA,xb,BAa12srscasssgg)

      call hp_a12carscrsssggg(xa,xb,ABa12carscrsssggg)
      call hp_a12carscrsssggg(xb,xa,BAa12carscrsssggg)
      call hp_a12carsssggg(xa,xb,ABa12carsssggg)
      call hp_a12carsssggg(xb,xa,BAa12carsssggg)
      call hp_a12carscasssggg(xa,xb,ABa12carscasssggg)
      call hp_a12carscasssggg(xb,xa,BAa12carscasssggg)
            

      tmpab3=
     1     half*ABa12srscrsgg-ABa12carssrscrsgg-two*ABa12carssrsssgg
     2     +ABa12srsssgg+two*ABa12carssrscasssgg-two*ABa12srscasssgg
     3     +two*ABa12carsssggg-two*ABa12carscasssggg-two*ABa12carscrsssggg
      tmpba3=
     1     half*BAa12srscrsgg-BAa12carssrscrsgg-two*BAa12carssrsssgg
     2     +BAa12srsssgg+two*BAa12carssrscasssgg-two*BAa12srscasssgg
     3     +two*BAa12carsssggg-two*BAa12carscasssggg-two*BAa12carscrsssggg
      
      tmpab3=-half*tmpab3
      tmpba3=-half*tmpba3
      
      k=0
      call hp_addscale2(k,tmp3,tmpab3,xl12)
      call hp_fillIop2_1(g,g,g,g,k,tmp3)

      call hp_addscale2(k,tmp3,tmpba3,xl12)
      call hp_fillIop2_2(g,g,g,g,k,tmp3)

c---  A1A1 Integrated counterterms
      
      call hp_a1a1carcasggg(xa,xb,ABa1a1carcasggg)
      call hp_a1a1carcasggg(xb,xa,BAa1a1carcasggg)
      call hp_a1a1cbrcasgggg(xa,xb,ABa1a1cbrcasgggg)
      call hp_a1a1cbrcasgggg(xb,xa,BAa1a1cbrcasgggg)
      call hp_a1a1crscasggg(xa,xb,ABa1a1crscasggg)
      call hp_a1a1crscasggg(xb,xa,BAa1a1crscasggg)

      tmpab2=
     1     two*ABa1a1carcasggg+two*ABa1a1cbrcasgggg+ABa1a1crscasggg

      tmpba2=
     1     two*BAa1a1carcasggg+two*BAa1a1cbrcasgggg+BAa1a1crscasggg   
 
      tmpab2=half*tmpab2
      tmpba2=half*tmpba2
           
      k=1
      call hp_addscale2(k,tmp2,tmpab2,xl12)
      call hp_fillIop2_1(g,g,g,g,k,tmp2)

      call hp_addscale2(k,tmp2,tmpba2,xl12)
      call hp_fillIop2_2(g,g,g,g,k,tmp2)
      
c---

      call hp_a1a1cassscarsragg(CA,xa,ABa1a1cassscarsragg)
      call hp_a1a1cassscarsragg(CA,xb,BAa1a1cassscarsragg)
      call hp_a1a1cassssragg(CA,xa,ABa1a1cassssragg)
      call hp_a1a1cassssragg(CA,xb,BAa1a1cassssragg)
      call hp_a1a1sscarsragg(CA,xa,ABa1a1sscarsragg)
      call hp_a1a1sscarsragg(CA,xb,BAa1a1sscarsragg)
      call hp_a1a1sssrgg(CA,xa,ABa1a1sssrgg)
      call hp_a1a1sssrgg(CA,xb,BAa1a1sssrgg)      

      call hp_a1a1crssscarggg(xa,xb,ABa1a1crssscarggg)
      call hp_a1a1crssscarggg(xb,xa,BAa1a1crssscarggg)
      call hp_a1a1sscarggg(xa,xb,ABa1a1sscarggg)
      call hp_a1a1sscarggg(xb,xa,BAa1a1sscarggg)

      tmpab3=two*ABa1a1sscarggg-six*ABa1a1crssscarggg    
     1     +ABa1a1sssrgg-two*ABa1a1sscarsragg
     2     -three*ABa1a1cassssragg+six*ABa1a1cassscarsragg  

      tmpba3=two*BAa1a1sscarggg-six*BAa1a1crssscarggg
     1     +BAa1a1sssrgg-two*BAa1a1sscarsragg
     2     -three*BAa1a1cassssragg+six*BAa1a1cassscarsragg
      
      tmpab3=half*tmpab3
      tmpba3=half*tmpba3

      k=0
      call hp_addscale2(k,tmp3,tmpab3,xl12)
      call hp_fillIop2_1(g,g,g,g,k,tmp3)

      call hp_addscale2(k,tmp3,tmpba3,xl12)
      call hp_fillIop2_2(g,g,g,g,k,tmp3)


c$$$      write(6,*) 'xa=',xa
c$$$      write(6,*) 'xb=',xb
c$$$      write(6,*)
c$$$c      write(6,*) 'I10'
c$$$      write(6,*)
c$$$      write(6,FMT3) 'I10opeps[1,1]={',I10op1eps(g,g,g,1,1,:)+I10op2eps(g,g,g,1,1,:)
c$$$      write(6,*)
c$$$      write(6,FMT3) 'I10opeps[1,2]={',I10op1eps(g,g,g,1,2,:)+I10op2eps(g,g,g,2,1,:)
c$$$      write(6,*)
c$$$      write(6,FMT3) 'I10opeps[2,1]={',I10op1eps(g,g,g,2,1,:)+I10op2eps(g,g,g,1,2,:)
c$$$      write(6,*)
c$$$      write(6,FMT3) 'I10opeps[2,2]={',I10op1eps(g,g,g,2,2,:)+I10op2eps(g,g,g,2,2,:)
c$$$      write(6,*)
c$$$      write(6,*)
c$$$      write(6,*) 'I20'
c$$$      write(6,*)
c$$$      write(6,FMT3) 'I20opeps[1,0]={',I20op1eps(g,g,g,g,1,0,:)+I20op2eps(g,g,g,g,0,1,:)
c$$$      write(6,*)                   
c$$$      write(6,FMT3) 'I20opeps[0,1]={',I20op1eps(g,g,g,g,0,1,:)+I20op2eps(g,g,g,g,1,0,:)
c$$$      write(6,*)                   
c$$$      write(6,FMT3) 'I20opeps[0,0]={',I20op1eps(g,g,g,g,0,0,:)+I20op2eps(g,g,g,g,0,0,:)
c$$$c      write(6,*) 'I20op1eps[0,0]=',I20op1eps(g,g,g,g,0,0,:)
c$$$c      write(6,*) 'I20op2eps[0,0]=',I20op2eps(g,g,g,g,0,0,:)
c$$$      write(6,*)                   
c$$$      write(6,FMT3) 'I20opeps[1,1]={',I20op1eps(g,g,g,g,1,1,:)+I20op2eps(g,g,g,g,1,1,:)
c$$$      write(6,*)                   
c$$$      write(6,FMT3) 'I20opeps[1,2]={',I20op1eps(g,g,g,g,1,2,:)+I20op2eps(g,g,g,g,2,1,:)
c$$$      write(6,*)
c$$$      write(6,FMT3) 'I20opeps[2,1]={',I20op1eps(g,g,g,g,2,1,:)+I20op2eps(g,g,g,g,1,2,:)
c$$$      write(6,*)
c$$$      write(6,FMT3) 'I20opeps[2,2]={',I20op1eps(g,g,g,g,2,2,:)+I20op2eps(g,g,g,g,2,2,:)
c$$$      write(6,*) 'in Iopmat'
c$$$      pause



      
c---  RVA1 Integrated counterterms

      call hp_rva1ncar1gg(xa,xb,ABrva1ncar1gg)
      call hp_rva1ncar1gg(xb,xa,BArva1ncar1gg)

      tmpab2=ABrva1ncar1gg

      tmpba2=BArva1ncar1gg

      k=1
      call hp_addscale2(k,tmp2,tmpab2,xl12)
      call hp_fillIop2_1(g,g,g,g,k,tmp2)

      call hp_addscale2(k,tmp2,tmpba2,xl12)
      call hp_fillIop2_2(g,g,g,g,k,tmp2)
      
c---
      
      call hp_rva1sr1g(CA,xa,ABrva1sr1g)
      call hp_rva1sr1g(CA,xb,BArva1sr1g)
      call hp_rva1carsr1g(CA,xa,ABrva1carsr1g)
      call hp_rva1carsr1g(CA,xb,BArva1carsr1g)

      tmpab3=half*ABrva1sr1g-ABrva1carsr1g

      tmpba3=half*BArva1sr1g-BArva1carsr1g

      k=0
      call hp_addscale2(k,tmp3,tmpab3,xl12)
      call hp_fillIop2_1(g,g,g,g,k,tmp3)

      call hp_addscale2(k,tmp3,tmpba3,xl12)
      call hp_fillIop2_2(g,g,g,g,k,tmp3)
      
      
c---  A1 Integrated counterterms
      
      call hp_a1ncar0gg(xa,xb,ABa1ncar0gg)
      call hp_a1ncar0gg(xb,xa,BAa1ncar0gg)


c---  !!! Note the following terms cancels, but could not be included without splitting the
c---  !!! lam,lam and collecting them into the 0,0 terms    
c      call hp_a1sr0g(CA,xa,ABa1sr0g)
c      call hp_a1sr0g(CA,xb,BAa1sr0g)

c      call hp_a1carsr0g(CA,xa,ABa1carsr0g)
c      call hp_a1carsr0g(CA,xb,BAa1carsr0g)

      tmpab1=ABa1ncar0gg!+half*ABa1sr0g-ABa1carsr0g

      tmpba1=BAa1ncar0gg!+half*BAa1sr0g-BAa1carsr0g
      
      call hp_addscale1(tmp1,tmpab1,xl12)
      call hp_fillIop1_1(g,g,g,1,tmp1)

      call hp_addscale1(tmp1,tmpba1,xl12)
      call hp_fillIop1_2(g,g,g,1,tmp1)
      


c      write(6,*) 'xa=',xa
c      write(6,*) 'xb=',xb
c      write(6,*)
c      write(6,*)
c      write(6,FMT3) 'I10op1eps[1,1]={',I10op1eps(g,g,g,1,1,:)+I10op2eps(g,g,g,1,1,:)
c      write(6,*)
c      write(6,FMT3) 'I10op1eps[1,2]={',I10op1eps(g,g,g,1,2,:)+I10op2eps(g,g,g,2,1,:)
c      write(6,*)
c      write(6,FMT3) 'I10op1eps[2,1]={',I10op1eps(g,g,g,2,1,:)+I10op2eps(g,g,g,1,2,:)
c      write(6,*)
c      write(6,FMT3) 'I10op1eps[2,2]={',I10op1eps(g,g,g,2,2,:)+I10op2eps(g,g,g,2,2,:)
c      write(6,*)
c      write(6,*)
c      pause



c$$$      write(6,*) 'xa=',xa
c$$$      write(6,*) 'xb=',xb
c$$$      write(6,*)
c$$$c      write(6,*) 'I10'
c$$$      write(6,*)
c$$$      write(6,FMT3) 'I10opeps[1,1]={',I10op1eps(g,g,g,1,1,:)+I10op2eps(g,g,g,1,1,:)
c$$$      write(6,*)
c$$$      write(6,FMT3) 'I10opeps[1,2]={',I10op1eps(g,g,g,1,2,:)+I10op2eps(g,g,g,2,1,:)
c$$$      write(6,*)
c$$$      write(6,FMT3) 'I10opeps[2,1]={',I10op1eps(g,g,g,2,1,:)+I10op2eps(g,g,g,1,2,:)
c$$$      write(6,*)
c$$$      write(6,FMT3) 'I10opeps[2,2]={',I10op1eps(g,g,g,2,2,:)+I10op2eps(g,g,g,2,2,:)
c$$$      write(6,*)
c$$$      write(6,*)
c$$$      write(6,*) 'I20'
c$$$      write(6,*)
c$$$      write(6,FMT3) 'I20opeps[1,0]={',I20op1eps(g,g,g,g,1,0,:)+I20op2eps(g,g,g,g,0,1,:)
c$$$      write(6,*)                   
c$$$      write(6,FMT3) 'I20opeps[0,1]={',I20op1eps(g,g,g,g,0,1,:)+I20op2eps(g,g,g,g,1,0,:)
c$$$      write(6,*)                   
c$$$      write(6,FMT3) 'I20opeps[0,0]={',I20op1eps(g,g,g,g,0,0,:)+I20op2eps(g,g,g,g,0,0,:)
c$$$c      write(6,*) 'I20op1eps[0,0]=',I20op1eps(g,g,g,g,0,0,:)
c$$$c      write(6,*) 'I20op2eps[0,0]=',I20op2eps(g,g,g,g,0,0,:)
c$$$      write(6,*)                   
c$$$      write(6,FMT3) 'I20opeps[1,1]={',I20op1eps(g,g,g,g,1,1,:)+I20op2eps(g,g,g,g,1,1,:)
c$$$      write(6,*)                   
c$$$      write(6,FMT3) 'I20opeps[1,2]={',I20op1eps(g,g,g,g,1,2,:)+I20op2eps(g,g,g,g,2,1,:)
c$$$      write(6,*)
c$$$      write(6,FMT3) 'I20opeps[2,1]={',I20op1eps(g,g,g,g,2,1,:)+I20op2eps(g,g,g,g,1,2,:)
c$$$      write(6,*)
c$$$      write(6,FMT3) 'I20opeps[2,2]={',I20op1eps(g,g,g,g,2,2,:)+I20op2eps(g,g,g,g,2,2,:)
c$$$      write(6,*) 'in Iopmat'
c$$$      pause



      
c--- Renormalization of the I2

      tmpab1=-tmpab1*b0
      tmpba1=-tmpba1*b0
      
      call hp_addscale3(tmp2,tmpab1,-xl12)
      call hp_fillIop2_1(g,g,g,g,1,tmp2)

      call hp_addscale3(tmp2,tmpba1,-xl12)
      call hp_fillIop2_2(g,g,g,g,1,tmp2)


c$$$      write(6,*) 'xa=',xa
c$$$      write(6,*) 'xb=',xb
c$$$      write(6,*)
c$$$c      write(6,*) 'I10'
c$$$      write(6,*)
c$$$      write(6,FMT3) 'I10opeps[1,1]={',I10op1eps(g,g,g,1,1,:)+I10op2eps(g,g,g,1,1,:)
c$$$      write(6,*)
c$$$      write(6,FMT3) 'I10opeps[1,2]={',I10op1eps(g,g,g,1,2,:)+I10op2eps(g,g,g,2,1,:)
c$$$      write(6,*)
c$$$      write(6,FMT3) 'I10opeps[2,1]={',I10op1eps(g,g,g,2,1,:)+I10op2eps(g,g,g,1,2,:)
c$$$      write(6,*)
c$$$      write(6,FMT3) 'I10opeps[2,2]={',I10op1eps(g,g,g,2,2,:)+I10op2eps(g,g,g,2,2,:)
c$$$      write(6,*)
c$$$      write(6,*)
c$$$      write(6,*) 'I20'
c$$$      write(6,*)
c$$$      write(6,FMT3) 'I20opeps[1,0]={',I20op1eps(g,g,g,g,1,0,:)+I20op2eps(g,g,g,g,0,1,:)
c$$$      write(6,*)                   
c$$$      write(6,FMT3) 'I20opeps[0,1]={',I20op1eps(g,g,g,g,0,1,:)+I20op2eps(g,g,g,g,1,0,:)
c$$$      write(6,*)                   
c$$$      write(6,FMT3) 'I20opeps[0,0]={',I20op1eps(g,g,g,g,0,0,:)+I20op2eps(g,g,g,g,0,0,:)
c$$$c      write(6,*) 'I20op1eps[0,0]=',I20op1eps(g,g,g,g,0,0,:)
c$$$c      write(6,*) 'I20op2eps[0,0]=',I20op2eps(g,g,g,g,0,0,:)
c$$$      write(6,*)                   
c$$$      write(6,FMT3) 'I20opeps[1,1]={',I20op1eps(g,g,g,g,1,1,:)+I20op2eps(g,g,g,g,1,1,:)
c$$$      write(6,*)                   
c$$$      write(6,FMT3) 'I20opeps[1,2]={',I20op1eps(g,g,g,g,1,2,:)+I20op2eps(g,g,g,g,2,1,:)
c$$$      write(6,*)
c$$$      write(6,FMT3) 'I20opeps[2,1]={',I20op1eps(g,g,g,g,2,1,:)+I20op2eps(g,g,g,g,1,2,:)
c$$$      write(6,*)
c$$$      write(6,FMT3) 'I20opeps[2,2]={',I20op1eps(g,g,g,g,2,2,:)+I20op2eps(g,g,g,g,2,2,:)
c$$$      write(6,*) 'in Iopmat'
c$$$      pause

      
c---  C1A1 Integrated counterterms
      
      call hp_a1c1cargggg(xa,xb,ABa1c1cargggg)
      call hp_a1c1cargggg(xb,xa,BAa1c1cargggg)
      call hp_a1c1cbrgggg(xa,xb,ABa1c1cbrgggg)
      call hp_a1c1cbrgggg(xb,xa,BAa1c1cbrgggg)

      tmpab1=ABa1c1cargggg+ABa1c1cbrgggg
      tmpba1=BAa1c1cargggg+BAa1c1cbrgggg

c---  xlfs=Log[Q^2/muF^2]
      
      call hp_addscale3(tmp,tmpab1,xlfs)
      call hp_addscale2(1,tmp2,tmp,xl12)
      call hp_fillIop2_1(g,g,g,g,1,tmp2)

      call hp_addscale3(tmp,tmpba1,xlfs)
      call hp_addscale2(1,tmp2,tmp,xl12)
      call hp_fillIop2_2(g,g,g,g,1,tmp2)


c$$$      write(6,*) 'xa=',xa
c$$$      write(6,*) 'xb=',xb
c$$$      write(6,*)
c$$$c      write(6,*) 'I10'
c$$$      write(6,*)
c$$$      write(6,FMT3) 'I10opeps[1,1]={',I10op1eps(g,g,g,1,1,:)+I10op2eps(g,g,g,1,1,:)
c$$$      write(6,*)
c$$$      write(6,FMT3) 'I10opeps[1,2]={',I10op1eps(g,g,g,1,2,:)+I10op2eps(g,g,g,2,1,:)
c$$$      write(6,*)
c$$$      write(6,FMT3) 'I10opeps[2,1]={',I10op1eps(g,g,g,2,1,:)+I10op2eps(g,g,g,1,2,:)
c$$$      write(6,*)
c$$$      write(6,FMT3) 'I10opeps[2,2]={',I10op1eps(g,g,g,2,2,:)+I10op2eps(g,g,g,2,2,:)
c$$$      write(6,*)
c$$$      write(6,*)
c$$$      write(6,*) 'I20'
c$$$      write(6,*)
c$$$      write(6,FMT3) 'I20opeps[1,0]={',I20op1eps(g,g,g,g,1,0,:)+I20op2eps(g,g,g,g,0,1,:)
c$$$      write(6,*)                   
c$$$      write(6,FMT3) 'I20opeps[0,1]={',I20op1eps(g,g,g,g,0,1,:)+I20op2eps(g,g,g,g,1,0,:)
c$$$      write(6,*)                   
c$$$      write(6,FMT3) 'I20opeps[0,0]={',I20op1eps(g,g,g,g,0,0,:)+I20op2eps(g,g,g,g,0,0,:)
c$$$c      write(6,*) 'I20op1eps[0,0]=',I20op1eps(g,g,g,g,0,0,:)
c$$$c      write(6,*) 'I20op2eps[0,0]=',I20op2eps(g,g,g,g,0,0,:)
c$$$      write(6,*)                   
c$$$      write(6,FMT3) 'I20opeps[1,1]={',I20op1eps(g,g,g,g,1,1,:)+I20op2eps(g,g,g,g,1,1,:)
c$$$      write(6,*)                   
c$$$      write(6,FMT3) 'I20opeps[1,2]={',I20op1eps(g,g,g,g,1,2,:)+I20op2eps(g,g,g,g,2,1,:)
c$$$      write(6,*)
c$$$      write(6,FMT3) 'I20opeps[2,1]={',I20op1eps(g,g,g,g,2,1,:)+I20op2eps(g,g,g,g,1,2,:)
c$$$      write(6,*)
c$$$      write(6,FMT3) 'I20opeps[2,2]={',I20op1eps(g,g,g,g,2,2,:)+I20op2eps(g,g,g,g,2,2,:)
c$$$      write(6,*) 'in Iopmat'
c$$$      pause

      
c--- AP C2 terms      

c---  xlrf=Log[muR^2/muF^2]

c---  First term of sigmaC2
c     asr/(2 Pi/2)*1/ep*(muR2/muF2)^(ep) (Pgg[0, 1 - xa]*OneLoop[xa pa, pb] DD[1 - xb] +  Pgg[0, 1 - xb]*OneLoop[pa, xb pb] DD[1 - xa])
            
      call hp_ap0(xa)
      ABPgg0=Pgg0
      call hp_ap0(xb)
      BAPgg0=Pgg0
      
c     Convert input function hp_with no eps dependence into eps-vector also adding
c     scale dependence (muR2/muF2)^(ep) and the 1/eps pole
      call hp_addscale1c(tmp1,ABPgg0,xlrf,1,1)
c     Accumulate the I10 operator which will be multiplied by the OneLoop ME
      call hp_fillIop1_1(g,g,g,1,tmp1)

c     Convert input function hp_with no eps dependence into eps-vector also adding
c     scale dependence (muR2/muF2)^(ep) and the 1/eps pole
      call hp_addscale1c(tmp1,BAPgg0,xlrf,1,1)
c     Accumulate the I10 operator which will be multiplied by the OneLoop ME
      call hp_fillIop1_2(g,g,g,1,tmp1)
            
c---  Second term of sigmaC2
c     (asr/(2 Pi/2))^2*1/(2 ep)*(muR2/muF2)^(2 ep) (Pgg[1, 1 - xa]*Born[xa pa, pb] DD[1 - xb] + Pgg[1, 1 - xb]*Born[pa, xb pb] DD[1 - xa])
           
      call hp_ap1(xa)
      ABPgg1=Pgg1
      call hp_ap1(xb)
      BAPgg1=Pgg1
            
c     Convert input function hp_with no eps dependence into eps-vector also adding
c     scale dependence (muR2/muF2)^(2 ep) and the 1/(2 eps) pole
      tempfcn = ABPgg1/2.0_ki
      call hp_addscale2c(tmp2,tempfcn,xlrf,1,2)
c     Accumulate the I20 operator which will be multiplied by the Born ME
      call hp_fillIop2_1(g,g,g,g,1,tmp2)

c     Convert input function hp_with no eps dependence into eps-vector also adding
c     scale dependence (muR2/muF2)^(2 ep) and the 1/(2 eps) pole
      tempfcn = BAPgg1/2.0_ki
      call hp_addscale2c(tmp2,tempfcn,xlrf,1,2)
c     Accumulate the I20 operator which will be multiplied by the Born ME
      call hp_fillIop2_2(g,g,g,g,1,tmp2)
      
c---  Third term of sigmaC2: part 1
c     (asr/(2 Pi/2))^2*beta0/(4 ep^2) ((muR2/muF2)^(2 ep) - 2 (muR2/muF2)^(ep)) (Pgg[0, 1 - xa]*Born[xa pa, pb] DD[1 - xb] + Pgg[0, 1 - xb]*Born[pa, xb pb] DD[1 - xa])      

c     Convert input function hp_with no eps dependence into eps-vector also adding
c     scale dependence (muR2/muF2)^(2 ep) and the beta0/(2 eps^2) pole
      tempfcn = ABPgg0*b0/2.0_ki
      call hp_addscale2c(tmp2,tempfcn,xlrf,2,2)
c     Accumulate the I20 operator which will be multiplied by the Born ME
      call hp_fillIop2_1(g,g,g,g,1,tmp2)

c     Convert input function hp_with no eps dependence into eps-vector also adding
c     scale dependence (muR2/muF2)^(2 ep) and the beta0/(2 eps^2) pole
      tempfcn = BAPgg0*b0/2.0_ki
      call hp_addscale2c(tmp2,tempfcn,xlrf,2,2)
c     Accumulate the I20 operator which will be multiplied by the Born ME
      call hp_fillIop2_2(g,g,g,g,1,tmp2)

c---  Third term of sigmaC2: part 2
c     (asr/(2 Pi/2))^2*beta0/(4 ep^2) ((muR2/muF2)^(2 ep) - 2 (muR2/muF2)^(ep)) (Pgg[0, 1 - xa]*Born[xa pa, pb] DD[1 - xb] + Pgg[0, 1 - xb]*Born[pa, xb pb] DD[1 - xa])      

c     Convert input function hp_with no eps dependence into eps-vector also adding
c     scale dependence (muR2/muF2)^(ep) and the -beta0/(eps^2) pole
      tempfcn = -ABPgg0*b0!/2.0_ki
      call hp_addscale2c(tmp2,tempfcn,xlrf,2,1)
c     Accumulate the I20 operator which will be multiplied by the Born ME
      call hp_fillIop2_1(g,g,g,g,1,tmp2)

c     Convert input function hp_with no eps dependence into eps-vector also adding
c     scale dependence (muR2/muF2)^(ep) and the -beta0/(eps^2) pole
      tempfcn = -BAPgg0*b0!/2.0_ki
      call hp_addscale2c(tmp2,tempfcn,xlrf,2,1)
c     Accumulate the I20 operator which will be multiplied by the Born ME
      call hp_fillIop2_2(g,g,g,g,1,tmp2)
      
c---  Fourth term of sigmaC2
c    (asr/(2 Pi/2))^2*1/(2 ep^2)*(muR2/muF2)^(2 ep) (ConvPxP[{g, g, 0}, {g, g, 0}, 1 - xa]*Born[xa pa, pb] DD[1 - xb] + ConvPxP[{g, g, 0}, {g, g, 0}, 1 - xb]*Born[pa, xb pb] DD[1 - xa])
      
      call hp_apxap(xa)
      ABP0ggxP0gg=P0ggxP0gg
      call hp_apxap(xb)
      BAP0ggxP0gg=P0ggxP0gg

c     Convert input function hp_with no eps dependence into eps-vector also adding
c     scale dependence (muR2/muF2)^(2 ep) and the 1/(2 eps^2) pole
      tempfcn = ABP0ggxP0gg/2.0_ki
      call hp_addscale2c(tmp2,tempfcn,xlrf,2,2)
c     Accumulate the I20 operator which will be multiplied by the Born ME
      call hp_fillIop2_1(g,g,g,g,1,tmp2)

c     Convert input function hp_with no eps dependence into eps-vector also adding
c     scale dependence (muR2/muF2)^(2 ep) and the 1/(2 eps^2) pole
      tempfcn = BAP0ggxP0gg/2.0_ki
      call hp_addscale2c(tmp2,tempfcn,xlrf,2,2)
c     Accumulate the I20 operator which will be multiplied by the Born ME
      call hp_fillIop2_2(g,g,g,g,1,tmp2)
      
c---  Fifth term of sigmaC2
c     (asr/(2 Pi/2))^2*1/ep^2*(muR2/muF2)^(2 ep) Pgg[0, 1 - xa]*Pgg[0, 1 - xb]*Born[xa pa, xb pb]

      call hp_apap(xa,xb)
      ABP0ggP0gg=P0ggP0gg
      call hp_apap(xb,xa)
      BAP0ggP0gg=P0ggP0gg

c     Convert input function hp_with no eps dependence into eps-vector also adding
c     scale dependence (muR2/muF2)^(2 ep) and the 1/(eps^2) pole
      tempfcn = ABP0ggP0gg/2.0_ki
      call hp_addscale2c(tmp2,tempfcn,xlrf,2,2)
c     Accumulate the I20 operator which will be multiplied by the Born ME
      call hp_fillIop2_1(g,g,g,g,1,tmp2)

c     Convert input function hp_with no eps dependence into eps-vector also adding
c     scale dependence (muR2/muF2)^(2 ep) and the 1/(eps^2) pole
      tempfcn = BAP0ggP0gg/2.0_ki
      call hp_addscale2c(tmp2,tempfcn,xlrf,2,2)
c     Accumulate the I20 operator which will be multiplied by the Born ME
      call hp_fillIop2_2(g,g,g,g,1,tmp2)



      
c$$$      write(6,*) 'xa=',xa
c$$$      write(6,*) 'xb=',xb
c$$$      write(6,*)
c$$$      write(6,*) 'I10'
c$$$      write(6,*)
c$$$      write(6,FMT3) 'I10opeps[1,1]={',I10op1eps(g,g,g,1,1,:)+I10op2eps(g,g,g,1,1,:)
c$$$      write(6,*)
c$$$      write(6,FMT3) 'I10opeps[1,2]={',I10op1eps(g,g,g,1,2,:)+I10op2eps(g,g,g,2,1,:)
c$$$      write(6,*)
c$$$      write(6,FMT3) 'I10opeps[2,1]={',I10op1eps(g,g,g,2,1,:)+I10op2eps(g,g,g,1,2,:)
c$$$      write(6,*)
c$$$      write(6,FMT3) 'I10opeps[2,2]={',I10op1eps(g,g,g,2,2,:)+I10op2eps(g,g,g,2,2,:)
c$$$      write(6,*)
c$$$      write(6,*)
c$$$      write(6,*) 'I20'
c$$$      write(6,*)
c$$$      write(6,FMT3) 'I20opeps[1,0]={',I20op1eps(g,g,g,g,1,0,:)+I20op2eps(g,g,g,g,0,1,:)
c$$$      write(6,*)                   
c$$$      write(6,FMT3) 'I20opeps[0,1]={',I20op1eps(g,g,g,g,0,1,:)+I20op2eps(g,g,g,g,1,0,:)
c$$$      write(6,*)                   
c$$$      write(6,FMT3) 'I20opeps[0,0]={',I20op1eps(g,g,g,g,0,0,:)+I20op2eps(g,g,g,g,0,0,:)
c$$$      write(6,*) 'I20op1eps[0,0]=',I20op1eps(g,g,g,g,0,0,:)
c$$$      write(6,*) 'I20op2eps[0,0]=',I20op2eps(g,g,g,g,0,0,:)
c$$$      write(6,*)                   
c$$$      write(6,FMT3) 'I20opeps[1,1]={',I20op1eps(g,g,g,g,1,1,:)+I20op2eps(g,g,g,g,1,1,:)
c$$$      write(6,*)                   
c$$$      write(6,FMT3) 'I20opeps[1,2]={',I20op1eps(g,g,g,g,1,2,:)+I20op2eps(g,g,g,g,2,1,:)
c$$$      write(6,*)
c$$$      write(6,FMT3) 'I20opeps[2,1]={',I20op1eps(g,g,g,g,2,1,:)+I20op2eps(g,g,g,g,1,2,:)
c$$$      write(6,*)
c$$$      write(6,FMT3) 'I20opeps[2,2]={',I20op1eps(g,g,g,g,2,2,:)+I20op2eps(g,g,g,g,2,2,:)
c$$$      pause

c 20   continue
      
      I20op1eps(:,:,:,:,:,:,:)=I20op1eps(:,:,:,:,:,:,:)!*ason2piprec2
      I20op2eps(:,:,:,:,:,:,:)=I20op2eps(:,:,:,:,:,:,:)!*ason2piprec2

      I10op1eps(:,:,:,:,:,:)=I10op1eps(:,:,:,:,:,:)!*ason2piprec
      I10op2eps(:,:,:,:,:,:)=I10op2eps(:,:,:,:,:,:)!*ason2piprec

      
      return
      end

      subroutine hp_fillIop2_1(a,b,c,d,i,Iopin)
      implicit none
      include 'hp_types.h'
      include 'hp_Ioperators2.f'
c      real(ki) 
c     1     I20op1eps(-1:1,-1:1,-1:1,-1:1,0:2,0:2,-4:0),
c     2     I20op2eps(-1:1,-1:1,-1:1,-1:1,0:2,0:2,-4:0)
c      common/I20opeps/I20op1eps,I20op2eps
      integer a,b,c,d,i,j,k,eps
      real(ki) Iopin(i:2,i:2,-4:0)
      do j=i,2
         do k=i,2
            do eps=-4,0
               I20op1eps(a,b,c,d,j,k,eps)=
     1              I20op1eps(a,b,c,d,j,k,eps)
     2              +Iopin(j,k,eps)
            enddo
         enddo
      enddo
      return
      end
      
      subroutine hp_fillIop2_2(a,b,c,d,i,Iopin)
      implicit none
      include 'hp_types.h'
      include 'hp_Ioperators2.f'
c      real(ki) 
c     1     I20op1eps(-1:1,-1:1,-1:1,-1:1,0:2,0:2,-4:0),
c     2     I20op2eps(-1:1,-1:1,-1:1,-1:1,0:2,0:2,-4:0)
c      common/I20opeps/I20op1eps,I20op2eps
      integer a,b,c,d,i,j,k,eps
      real(ki) Iopin(i:2,i:2,-4:0)
      do j=i,2
         do k=i,2
            do eps=-4,0
               I20op2eps(a,b,c,d,j,k,eps)=
     1              I20op2eps(a,b,c,d,j,k,eps)
     2              +Iopin(j,k,eps)
            enddo
         enddo
      enddo
      return
      end
      
      subroutine hp_fillIop1_1(a,b,c,i,Iopin)
      implicit none
      include 'hp_types.h'
      include 'hp_Ioperators2.f'
c      real(ki) 
c     1     I10op1eps(-1:1,-1:1,-1:1,2,2,-2:2),
c     2     I10op2eps(-1:1,-1:1,-1:1,2,2,-2:2)
c      common/I10opeps/I10op1eps,I10op2eps
      integer a,b,c,i,j,k,eps
      real(ki) Iopin(i:2,i:2,-2:2)
      do j=i,2
         do k=i,2
            do eps=-2,2
               I10op1eps(a,b,c,j,k,eps)=
     1              I10op1eps(a,b,c,j,k,eps)
     2              +Iopin(j,k,eps)
            enddo
         enddo
      enddo
      return
      end
      
      subroutine hp_fillIop1_2(a,b,c,i,Iopin)
      implicit none
      include 'hp_types.h'
      include 'hp_Ioperators2.f'
c      real(ki) 
c     1     I10op1eps(-1:1,-1:1,-1:1,2,2,-2:2),
c     2     I10op2eps(-1:1,-1:1,-1:1,2,2,-2:2)
c      common/I10opeps/I10op1eps,I10op2eps
      integer a,b,c,i,j,k,eps
      real(ki) Iopin(i:2,i:2,-2:2)
      do j=i,2
         do k=i,2
            do eps=-2,2
               I10op2eps(a,b,c,j,k,eps)=
     1              I10op2eps(a,b,c,j,k,eps)
     2              +Iopin(j,k,eps)
            enddo
         enddo
      enddo
      return
      end

      subroutine hp_addscale3(outarr,inarr,xlfs)
      implicit none
      include 'hp_types.h'
      include 'hp_constants.h'
      integer i,j
      real(ki) outarr(2,2,-4:0),inarr(2,2,-2:2)
      real(ki) xlfs
      do i=1,2
         do j=1,2
            outarr(i,j,-4)=0.0_ki
            outarr(i,j,-3)=inarr(i,j,-2)
            outarr(i,j,-2)=inarr(i,j,-1)
     1           +inarr(i,j,-2)*xlfs
            outarr(i,j,-1)=inarr(i,j, 0)
     1           +inarr(i,j,-1)*xlfs
     2           +inarr(i,j,-2)*(xlfs**2+pisq/6.0_ki)/2.0_ki
            outarr(i,j, 0)=inarr(i,j, 1)
     1           +inarr(i,j, 0)*xlfs
     2           +inarr(i,j,-1)*(xlfs**2+pisq/6.0_ki)/2.0_ki
     3           +inarr(i,j,-2)*(xlfs*pisq/6.0_ki+xlfs**3/3.0_ki+2*zeta3/3.0_ki)/2.0_ki
         enddo
      enddo
      return
      end

      subroutine hp_addscale2(k,outarr,inarr,xl12)
      implicit none
      include 'hp_types.h'
      integer i,j,k
      real(ki) inarr(k:2,k:2,-4:0),outarr(k:2,k:2,-4:0)
      real(ki) xl12
      do i=k,2
         do j=k,2
            outarr(i,j,-4)=inarr(i,j,-4)
            outarr(i,j,-3)=inarr(i,j,-3)
     1           -inarr(i,j,-4)*2*xl12
            outarr(i,j,-2)=inarr(i,j,-2)
     1           -inarr(i,j,-3)*2*xl12
     2           +inarr(i,j,-4)*2*xl12**2
            outarr(i,j,-1)=inarr(i,j,-1)
     1           -inarr(i,j,-2)*2*xl12
     2           +inarr(i,j,-3)*2*xl12**2
     3           -inarr(i,j,-4)*4/3.0_ki*xl12**3
            outarr(i,j, 0)=inarr(i,j, 0)
     1           -inarr(i,j,-1)*2*xl12
     2           +inarr(i,j,-2)*2*xl12**2
     3           -inarr(i,j,-3)*4/3.0_ki*xl12**3
     4           +inarr(i,j,-4)*2/3.0_ki*xl12**4
         enddo
      enddo
      return
      end
      
      subroutine hp_addscale1(outarr,inarr,xl12)
      implicit none
      include 'hp_types.h'
      integer i,j
      real(ki) inarr(2,2,-2:2),outarr(2,2,-2:2)
      real(ki) xl12
      do i=1,2
         do j=1,2
            outarr(i,j,-2)=inarr(i,j,-2)
            outarr(i,j,-1)=inarr(i,j,-1)
     1           -inarr(i,j,-2)*xl12
            outarr(i,j, 0)=inarr(i,j, 0)
     1           -inarr(i,j,-1)*xl12
     2           +inarr(i,j,-2)/2.0_ki*xl12**2
            outarr(i,j, 1)=inarr(i,j, 1)
     1           -inarr(i,j, 0)*xl12
     2           +inarr(i,j,-1)/2.0_ki*xl12**2
     3           -inarr(i,j,-2)/6.0_ki*xl12**3
            outarr(i,j, 2)=inarr(i,j, 2)
     1           -inarr(i,j, 1)*xl12
     2           +inarr(i,j, 0)/2.0_ki*xl12**2
     3           -inarr(i,j,-1)/6.0_ki*xl12**3
     4           +inarr(i,j,-2)/24.0_ki*xl12**4
         enddo
      enddo
      return
      end

      subroutine hp_addscale2c(outarr,infunc,xlrs,epspoleorder,logpower)
      implicit none
      include 'hp_types.h'
      include 'hp_constants.h'
      integer i,j, epspoleorder, logpower
      real(ki) infunc(2,2),outarr(2,2,-4:0)
      real(ki) xlrs
      if ((epspoleorder.ne.1).and.(epspoleorder.ne.2)) then
         write(6,*) "Unimplmented epsilon oder in addscale2c"
         stop
      endif
      do i=1,2
         do j=1,2
            if (epspoleorder.eq.1) then
               outarr(i,j,-4) = 0.0_ki
               outarr(i,j,-3) = 0.0_ki
               outarr(i,j,-2) = 0.0_ki
               outarr(i,j,-1) = infunc(i,j)
               outarr(i,j, 0) = infunc(i,j)*logpower*xlrs
            elseif (epspoleorder.eq.2) then
               outarr(i,j,-4) = 0.0_ki
               outarr(i,j,-3) = 0.0_ki
               outarr(i,j,-2) = infunc(i,j)
               outarr(i,j,-1) = infunc(i,j)*logpower*xlrs
               outarr(i,j, 0) = infunc(i,j)*((logpower*xlrs)**2/2.0_ki + pisq/6.0_ki)
            endif
         enddo
      enddo
      return
      end
      
      subroutine hp_addscale1c(outarr,infunc,xlrs,epspoleorder,logpower)
      implicit none
      include 'hp_types.h'
      include 'hp_constants.h'
      integer i,j, epspoleorder, logpower
      real(ki) infunc(2,2),outarr(2,2,-2:2)
      real(ki) xlrs
      if ((epspoleorder.ne.1).and.(epspoleorder.ne.2)) then
         write(6,*) "Unimplmented epsilon oder in addscale1c"
         stop
      endif
      do i=1,2
         do j=1,2
            if (epspoleorder.eq.1) then
               outarr(i,j,-2) = 0.0_ki
               outarr(i,j,-1) = infunc(i,j)
               outarr(i,j, 0) = infunc(i,j)*logpower*xlrs
               outarr(i,j, 1) = infunc(i,j)*((logpower*xlrs)**2/2.0_ki + pisq/12.0_ki)
               outarr(i,j, 2) = infunc(i,j)*((logpower*xlrs)**3/6.0_ki + pisq/12.0_ki*(logpower*xlrs) + zeta3/3.0_ki)
            elseif (epspoleorder.eq.2) then
               outarr(i,j,-2) = infunc(i,j)
               outarr(i,j,-1) = infunc(i,j)*logpower*xlrs
               outarr(i,j, 0) = infunc(i,j)*((logpower*xlrs)**2/2.0_ki + pisq/12.0_ki)
               outarr(i,j, 1) = infunc(i,j)*((logpower*xlrs)**3/6.0_ki + pisq/12.0_ki*(logpower*xlrs) + zeta3/3.0_ki)
               outarr(i,j, 2) = infunc(i,j)*((logpower*xlrs)**4/24.0_ki + pisq/24.0_ki*(logpower*xlrs)**2 + (logpower*xlrs)*zeta3/3.0_ki + pisq**2/160.0_ki)
            endif
         enddo
      enddo
      return
      end
      
      
