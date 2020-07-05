# SCCTB makefile
# SCCTB library rules
# Might modify optimization levels.
# The headfile is handled in install.com
# QC: change FC2 to FC to avoid free format (but need -c)
FFLAGS3 = -O3
FCSCC = $(FC) -c

OBJS_sccdftb= \
	atomen.o \
	broyden.o \
	cpe.o \
	dftd3.o \
	dispersionread.o \
	dispersion_egr.o \
	diis.o \
	dxlc.o \
	dylcao.o \
	eglcao.o \
	ewevge.o \
	externalchgrad.o \
	externalshift.o \
	fermi.o \
	gamma.o \
	gammamat.o \
	geometries.o \
	gettab.o \
	gradient.o \
	inicof.o \
	ktable_ewscc.o \
	long_range.o \
	long_sccmm1.o \
	long_sccmm2.o \
        long_sccpme.o \
	mixer.o \
	mulliken.o \
	output.o \
	repulsiv.o \
        sccpme.o \
	self.o \
	shift.o \
	short_range.o \
	skpar.o \
	slkode.o \
	slktrafo.o 
#
$(LIB)/sccdftb.a : $(OBJS_sccdftb)
	$(AR_COMMAND) $(LIB)/sccdftb.a $(OBJS_sccdftb)
	$(RANLIB) $(LIB)/sccdftb.a
	@echo SCCDFTB COMPLETED
#
# SCCTB source file rules
atomen.o : $(SRC)/sccdftbint/sccdftbsrc/atomen.f
	$(FCSCC) $(FFLAGS3) $(SRC)/sccdftbint/sccdftbsrc/atomen.f

broyden.o : $(SRC)/sccdftbint/sccdftbsrc/broyden.f
	$(FCSCC) $(FFLAGS3) $(SRC)/sccdftbint/sccdftbsrc/broyden.f

cpe.o : $(SRC)/sccdftbint/sccdftbsrc/cpe.F90
	$(FC2) $^

# sccdftbint source file rules
dftd3.o : $(SRC)/sccdftbint/sccdftbsrc/dftd3.F90
	$(FC2) $^

diis.o : $(SRC)/sccdftbint/sccdftbsrc/diis.F90
	$(FC2) $^

dispersionread.o : $(SRC)/sccdftbint/sccdftbsrc/dispersionread.f
	$(FCSCC) $(FFLAGS3) $(SRC)/sccdftbint/sccdftbsrc/dispersionread.f

dispersion_egr.o : $(SRC)/sccdftbint/sccdftbsrc/dispersion_egr.f
	$(FCSCC) $(FFLAGS3) $(SRC)/sccdftbint/sccdftbsrc/dispersion_egr.f

dxlc.o : $(SRC)/sccdftbint/sccdftbsrc/dxlc.f
	$(FCSCC) $(FFLAGS3) $(SRC)/sccdftbint/sccdftbsrc/dxlc.f

dylcao.o : $(SRC)/sccdftbint/sccdftbsrc/dylcao.F
	$(FCSCC) $(FFLAGS3) $^

eglcao.o : $(SRC)/sccdftbint/sccdftbsrc/eglcao.F
	$(FCSCC) $(FFLAGS3) $^

ewevge.o : $(SRC)/sccdftbint/sccdftbsrc/ewevge.F
	$(FCSCC) $(FFLAGS3) $^

externalchgrad.o : $(SRC)/sccdftbint/sccdftbsrc/externalchgrad.F
	$(FCSCC) $(FFLAGS3) $^

externalshift.o : $(SRC)/sccdftbint/sccdftbsrc/externalshift.f
	$(FCSCC) $(FFLAGS3) $(SRC)/sccdftbint/sccdftbsrc/externalshift.f

fermi.o : $(SRC)/sccdftbint/sccdftbsrc/fermi.f
	$(FCSCC) $(FFLAGS3) $(SRC)/sccdftbint/sccdftbsrc/fermi.f

gamma.o : $(SRC)/sccdftbint/sccdftbsrc/gamma.f
	$(FCSCC) $(FFLAGS3) $(SRC)/sccdftbint/sccdftbsrc/gamma.f

gammamat.o : $(SRC)/sccdftbint/sccdftbsrc/gammamat.f
	$(FCSCC) $(FFLAGS3) $(SRC)/sccdftbint/sccdftbsrc/gammamat.f

geometries.o : $(SRC)/sccdftbint/sccdftbsrc/geometries.f
	$(FCSCC) $(FFLAGS3) $(SRC)/sccdftbint/sccdftbsrc/geometries.f

gettab.o : $(SRC)/sccdftbint/sccdftbsrc/gettab.f
	$(FCSCC) $(FFLAGS3) $(SRC)/sccdftbint/sccdftbsrc/gettab.f

gradient.o : $(SRC)/sccdftbint/sccdftbsrc/gradient.F
	$(FCSCC) $(FFLAGS3) $^

inicof.o : $(SRC)/sccdftbint/sccdftbsrc/inicof.f
	$(FCSCC) $(FFLAGS3) $(SRC)/sccdftbint/sccdftbsrc/inicof.f

# Treat ktable with FLX 'cause we use the PARALLEL flag
ktable_ewscc.o : $(SRC)/sccdftbint/sccdftbsrc/ktable_ewscc.F
	$(FCSCC) $(FFLAGS3) $^

long_range.o : $(SRC)/sccdftbint/sccdftbsrc/long_range.f
	$(FCSCC) $(FFLAGS3) $(SRC)/sccdftbint/sccdftbsrc/long_range.f

long_sccmm1.o : $(SRC)/sccdftbint/sccdftbsrc/long_sccmm1.F
	$(FCSCC) $^

long_sccmm2.o : $(SRC)/sccdftbint/sccdftbsrc/long_sccmm2.F
	$(FCSCC) $(FFLAGS3) $^

long_sccpme.o : $(SRC)/sccdftbint/sccdftbsrc/long_sccpme.F
	$(FCSCC) $(FFLAGS3) $^

mixer.o : $(SRC)/sccdftbint/sccdftbsrc/mixer.F
	$(FCSCC) $(FFLAGS3) $^

mulliken.o : $(SRC)/sccdftbint/sccdftbsrc/mulliken.F
	$(FCSCC) $(FFLAGS3) $^

output.o : $(SRC)/sccdftbint/sccdftbsrc/output.f
	$(FCSCC) $(FFLAGS3) $(SRC)/sccdftbint/sccdftbsrc/output.f

repulsiv.o : $(SRC)/sccdftbint/sccdftbsrc/repulsiv.f
	$(FCSCC) $(FFLAGS3) $(SRC)/sccdftbint/sccdftbsrc/repulsiv.f

self.o : $(SRC)/sccdftbint/sccdftbsrc/self.f
	$(FCSCC) $(FFLAGS3) $(SRC)/sccdftbint/sccdftbsrc/self.f

shift.o : $(SRC)/sccdftbint/sccdftbsrc/shift.f
	$(FCSCC) $(FFLAGS3) $(SRC)/sccdftbint/sccdftbsrc/shift.f

short_range.o : $(SRC)/sccdftbint/sccdftbsrc/short_range.f
	$(FCSCC) $(FFLAGS3) $(SRC)/sccdftbint/sccdftbsrc/short_range.f

skpar.o : $(SRC)/sccdftbint/sccdftbsrc/skpar.f
	$(FCSCC) $(FFLAGS3) $(SRC)/sccdftbint/sccdftbsrc/skpar.f

slkode.o : $(SRC)/sccdftbint/sccdftbsrc/slkode.f
	$(FCSCC) $(FFLAGS3) $(SRC)/sccdftbint/sccdftbsrc/slkode.f

slktrafo.o : $(SRC)/sccdftbint/sccdftbsrc/slktrafo.f
	$(FCSCC) $(FFLAGS3) $(SRC)/sccdftbint/sccdftbsrc/slktrafo.f

#
# SCCDFTB dependency file
#

atomen.o : sccdftbsrc_ltm.o

broyden.o : sccdftbsrc_ltm.o

cpe.o : sccdftb_ltm.o sccdftbsrc_ltm.o psf_ltm.o

dftd3.o : sccdftb_ltm.o sccdftbsrc_ltm.o sccgsbp_ltm.o sccpb_ltm.o chm_kinds_ltm.o

dispersionread.o : sccdftbsrc_ltm.o

dispersion_egr.o : sccdftbsrc_ltm.o

dylcao.o : blockscc_ltm.o dimens_ltm.o coord_ltm.o sccdftb_ltm.o sccdftbsrc_ltm.o dftd3.o cpe.o

eglcao.o : sccdftb_ltm.o sccdftbsrc_ltm.o sccgsbp_ltm.o sccpb_ltm.o chm_kinds_ltm.o dftd3.o cpe.o diis.o  new_timer.o

externalchgrad.o : sccdftb_ltm.o sccdftbsrc_ltm.o

externalshift.o : sccdftb_ltm.o sccdftbsrc_ltm.o

fermi.o : sccdftbsrc_ltm.o

gammamat.o : erfcd.o sccdftb_ltm.o sccdftbsrc_ltm.o

geometries.o : sccdftb_ltm.o

gettab.o : sccdftbsrc_ltm.o

gradient.o : sccdftbsrc_ltm.o

ktable_ewscc.o : parallel_ltm.o

long_range.o : sccdftb_ltm.o sccdftbsrc_ltm.o

long_sccmm1.o : parallel_ltm.o sccdftb_ltm.o sccdftbsrc_ltm.o

long_sccmm2.o : parallel_ltm.o sccdftb_ltm.o sccdftbsrc_ltm.o

long_sccpme.o : consta_ltm.o exfunc_ltm.o grape_ltm.o memory_mod.o number_ltm.o sccpme.o sccpmeutil.o sccdftb_ltm.o sccdftbsrc_ltm.o 

mixer.o : sccdftbsrc_ltm.o

mulliken.o : sccdftb_ltm.o sccdftbsrc_ltm.o new_timer.o

output.o : sccdftbsrc_ltm.o

repulsiv.o : sccdftbsrc_ltm.o

shift.o : sccdftbsrc_ltm.o

short_range.o : sccdftbsrc_ltm.o

skpar.o : sccdftbsrc_ltm.o

slkode.o : sccdftb_ltm.o sccdftbsrc_ltm.o
