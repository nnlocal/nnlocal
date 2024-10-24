# Makefile routine.

# PATHS
NNLOCALHOME     = $(PWD)
SOURCEDIR       = $(NNLOCALHOME)/src
VPATH		= $(DIRS)
BIN		= $(NNLOCALHOME)/bin
INCPATH  	= $(SOURCEDIR)/Inc
INCPATH2  	= $(SOURCEDIR)/Intcts/include
INCPATH3  	= $(SOURCEDIR)/hp_Intcts/include
OUTPUT_OPTION	= -o $(NNLOCALHOME)/obj/$@

#LHAPDFLIB
PDFROUTINES = LHAPDF
LHAPDF_CONFIG=lhapdf-config
LHAPDFLIB = $(shell $(LHAPDF_CONFIG) --libdir)

#FASTJET_CONFIG=lhapdf-config
#FASTJETLIB = $(shell $(FASTJET_CONFIG) --libs)
#LIBFLAGS += $(FASTJETLIB)

FC = gfortran
FFLAGS 	= -fno-automatic -std=legacy -fno-f2c -g -ffixed-line-length-none -I$(INCPATH) -I$(INCPATH2) -I$(INCPATH3) -Iobj
#FFLAGS 	= -fno-automatic -std=legacy -fno-f2c -O2 -g -ffixed-line-length-none -I$(INCPATH) -I$(INCPATH2) -I$(INCPATH3) -Iobj

F90 = gfortran
F90FLAGS = -fno-automatic -fno-f2c  -g -I$(INCPATH) -I$(INCPATH2) -I$(INCPATH3) -Iobj -Jobj

CXX= c++
CXXFLAGS+= -std=c++17



DIRS	=	$(NNLOCALHOME):\
		$(NNLOCALHOME)/obj:\
		$(SOURCEDIR)/User:\
		$(SOURCEDIR)/Procdep:\
		$(SOURCEDIR)/Need:\
		$(SOURCEDIR)/Phase:\
		$(SOURCEDIR)/pslimits:\
		$(SOURCEDIR)/Intcts:\
		$(SOURCEDIR)/hp_Intcts:\
		$(SOURCEDIR)/Parton:\
		$(SOURCEDIR)/Integrate:\
                $(SOURCEDIR)/ggH:

# -----------------------------------------------------------------------------
# Specify the object files. 

HPINTCTS = \
hp_Iopmat2.o \
hp_cli2.o \
hp_Li3.o \
hp_WGPLG.o \
hp_mylog.o \
hp_glog.o \
hp_a12srscrsgg.o \
hp_a2srsgg.o \
hp_a2carssrsgg.o \
hp_a2carbssrsgg.o \
hp_a2carsggg.o \
hp_a2carsggg_bulk.o \
hp_a2carbsgggg.o \
hp_a2carbsgggg_bulk.o \
hp_a12srsssgg.o \
hp_a12srscasssgg.o \
hp_a12carscasggg.o \
hp_a12carssrsssgg.o \
hp_a12casbrcasgggg.o \
hp_a12carssrscrsgg.o \
hp_a12carssrscasssgg.o \
hp_AGTBa12carscrsggg.o \
hp_BGTAa12carscrsggg.o \
hp_a12carscrs.o \
hp_AGTBa12carscrsssggg.o \
hp_BGTAa12carscrsssggg.o \
hp_a12carscrsss.o \
hp_AGTBa12carscasssggg.o \
hp_BGTAa12carscasssggg.o \
hp_a12carscasss.o \
hp_AGTBa12carsssggg.o \
hp_BGTAa12carsssggg.o \
hp_a12carsss.o \
hp_a1a1cassscarsragg.o \
hp_a1a1cbrcasgggg.o \
hp_a1a1cassssragg.o \
hp_a1a1sscarsragg.o \
hp_a1a1carcasggg.o \
hp_a1a1sssrgg.o \
hp_AGTBa1a1crscasggg.o \
hp_BGTAa1a1crscasggg.o \
hp_a1a1crscas.o \
hp_AGTBa1a1crssscarggg.o \
hp_BGTAa1a1crssscarggg.o \
hp_a1a1crssscar.o \
hp_AGTBa1a1sscarggg.o \
hp_BGTAa1a1sscarggg.o \
hp_a1a1sscar.o \
hp_rva1ncar1gg.o \
hp_rva1sr1g.o \
hp_rva1carsr1g.o \
hp_a1ncar0gg.o \
hp_a1sr0g.o \
hp_a1carsr0g.o \
hp_a1c1cargggg.o \
hp_a1c1cbrgggg.o \
hp_APfactors.o


INTCTS = \
Iopmat2.o \
cli2.o \
Li3.o \
WGPLG.o \
mylog.o \
glog.o \
a12srscrsgg.o \
a2srsgg.o \
a2carssrsgg.o \
a2carbssrsgg.o \
a2carsggg.o \
a2carsggg_bulk.o \
a2carbsgggg.o \
a2carbsgggg_bulk.o \
a12srsssgg.o \
a12srscasssgg.o \
a12carscasggg.o \
a12carssrsssgg.o \
a12casbrcasgggg.o \
a12carssrscrsgg.o \
a12carssrscasssgg.o \
AGTBa12carscrsggg.o \
BGTAa12carscrsggg.o \
a12carscrs.o \
AGTBa12carscrsssggg.o \
BGTAa12carscrsssggg.o \
a12carscrsss.o \
AGTBa12carscasssggg.o \
BGTAa12carscasssggg.o \
a12carscasss.o \
AGTBa12carsssggg.o \
BGTAa12carsssggg.o \
a12carsss.o \
a1a1cassscarsragg.o \
a1a1cbrcasgggg.o \
a1a1cassssragg.o \
a1a1sscarsragg.o \
a1a1carcasggg.o \
a1a1sssrgg.o \
AGTBa1a1crscasggg.o \
BGTAa1a1crscasggg.o \
a1a1crscas.o \
AGTBa1a1crssscarggg.o \
BGTAa1a1crssscarggg.o \
a1a1crssscar.o \
AGTBa1a1sscarggg.o \
BGTAa1a1sscarggg.o \
a1a1sscar.o \
rva1ncar1gg.o \
rva1sr1g.o \
rva1carsr1g.o \
a1ncar0gg_loprec.o \
a1ncar0gg.o \
a1sr0g.o \
a1carsr0g.o \
a1c1cargggg.o \
a1c1cbrgggg.o \
APfactors.o


GGHFILES = \
finitemtcorr.o \
Ftriangle.o \
gg_hgaga_v_cf.o \
gg_hgamgam.o \
gg_hgamgam_v.o \
gg_hgamgam_gvec.o \
gg_hgamgam_gvecgvec.o \
gg_hgamgam_v_cf.o \
gg_hgamgam_gs_cf.o \
gg_hgamgam_vgvec_cf.o \
gg_hgagag_gs_colorful.o \
gg_hgagag.o \
gg_hgagagg.o \
gg_hgagag_gvec.o \
gg_hgagag_v.o \
gg_hg_z_cf.o \
gg_hg_z_wrap.o \
gg_htot.o \
h4g.o \
hqqgg.o \
hjetfill.o \
msqgamgam.o


INTEGRATEFILES = \
dgauss.o \
ebook.o \
mbook.o \
ran0.o \
ran1.o \
rn.o \
vegas.o


NEEDFILES = \
aveptjet.o \
banner.o \
boost.o \
checkndotp.o \
checkversion.o \
ckmfill.o \
compnewgridpar.o \
computepdfuncertainty.o \
coupling.o \
coupling2.o \
couplz.o \
ddilog.o \
dzero.o \
dipolesubA1_colorful.o \
dipolesubRVA1_colorful.o \
dipolesubA2_colorful.o \
dipolesubA12_colorful.o \
determinefilenames.o \
dipoles_cf.o \
donothing_gvec.o \
dot.o \
dotem.o \
etmiss.o \
getbs.o \
getptilde.o \
getptildejet.o \
higgsp.o \
higgsw.o \
histofin.o \
hbbdecay.o \
htautaudecay.o \
hwwdecay.o \
hbbdecay_g.o \
hbbdecay_v.o \
Iopmat1.o \
includedipole.o \
interpolate_hto.o \
is_functions.o \
lenocc.o \
lnrat.o \
masscuts.o \
mfrun.o \
nnlocal_main.o \
nnlocal_exit.o \
nnlocal_init.o \
nnlocal_vegas.o \
ptyrap.o \
r.o \
read_jetcuts.o \
reader_input.o \
realhistos.o \
scaleset.o \
scaleset_m34.o \
scaleset_Msqpt34sq.o \
scaleset_HT.o \
sethparams.o \
setmb_msbar.o \
setrunname.o \
setvdecay.o \
smalls.o \
spinoru.o \
storeptilde.o \
swapjet.o \
transformA1_colorful.o \
transformA2_colorful.o \
writeinfo.o \
writeinput.o \
writeout.o \
zeromsq.o


PARTONFILES = \
alfamz.o \
checkpath.o \
newton1.o


PSLIMITFILES = \
genlimit.o \
print.o \
rambo.o \
singen.o


PHASEFILES = \
breitw.o \
gencheck.o \
gen2.o \
gen3.o \
gen4.o \
phase3.o \
phase4.o \
phi1_2m_nobw.o \
phi3.o \
phi3m.o \
phi3m0.o


PROCDEPFILES = \
bornint.o \
chooser.o \
gen_born_ps.o \
gen_virt_ps.o \
gen_real_ps.o \
realint.o \
virtint.o


USERFILES = \
bookplot.o \
deltarj.o \
durhamalg.o \
genclust2.o \
genclust_kt.o \
genclust_cone.o \
gencuts.o \
gencuts_input.o \
genplots.o \
getet.o \
idjet.o \
integratehisto.o \
irregbins.o \
mdata.o \
miscclust.o \
nplotter.o \
nplotter_generic.o \
setnotag.o

LIBDIR=.

MAIN = nnlocal.o

ifeq ($(PDFROUTINES),LHAPDF)
   PARTONFILES += \
   alfamz_lhapdf.o \
   fdist_lhapdf.o \
   pdfwrap_lhapdf.o
   LIBDIR += -L$(LHAPDFLIB)
   LIBFLAGS += -lLHAPDF
   PDFMSG='   ----> NNLOCAL compiled with LHAPDF routines <----'
endif

OURCODE = $(NEEDFILES)  $(PROCDEPFILES) \
          $(PHASEFILES) $(PSLIMITFILES) $(HPINTCTS) $(INTCTS) \
          $(USERFILES) $(GGHFILES)

OTHER = $(INTEGRATEFILES) $(PARTONFILES)

ALLNNLOCAL = $(OURCODE) $(OTHER) $(MAIN)

obj:
	mkdir -p $@

nnlocal: $(ALLNNLOCAL)
	$(FC) $(FFLAGS) -L$(LIBDIR) -o $@ \
	$(patsubst %,obj/%,$(ALLNNLOCAL)) $(LIBFLAGS) 
	mv nnlocal bin/
	@echo $(PDFMSG)

gridplot: rgrids.o
	$(FC) $(FFLAGS) -o  $@ \
	$(patsubst %,obj/%,rgrids.o)
	mv gridplot bin/

mergedata: mergedata.o
	$(FC) $(FFLAGS) -o  $@ \
	$(patsubst %,obj/%,mergedata.o)
	mv mergedata bin/

%.o: %.f90
	$(F90) $(F90FLAGS) -c -o obj/$@ $<

%.o: %.cc
	$(CXX) -c $(CXXFLAGS) `fastjet-config --cxxflags` -o obj/$@ $<

# -----------------------------------------------------------------------------
# Specify other options.

# Specify the dependencies of the .o files and the rules to make them.

clean:
	- rm -f *.o obj/*.o obj/*.mod bin/nnlocal bin/gridplot bin/mergedata*.s *~ core

very-clean:
	- rm -f *.o bin/nnlocal bin/gridplot bin/mergedata *.s *~ core
	- rm -fr obj

# -----------------------------------------------------------------------------

all: obj nnlocal gridplot mergedata

.DEFAULT_GOAL:= all

# DO NOT DELETE

