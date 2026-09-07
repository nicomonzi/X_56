## $Header: /var/cvs/mbdyn/mbdyn/mbdyn-1.0/tests/rtai/pendulum-control/pendulum_control_rtai/pendulum_control.mk,v 1.3 2007/02/26 00:34:50 masarati Exp $

MAKECMD         = make
HOST=UNIX
BUILD=yes
SYS_TARGET_FILE = rtai.tlc

MODEL           = pendulum_control
MODULES         = pendulum_control_data.c 
MAKEFILE        = pendulum_control.mk
MATLAB_ROOT     = /usr/local/matlab
S_FUNCTIONS     = sfun_mbdyn_com_write.c sfun_mbdyn_com_read.c sfun_rtai_scope.c
S_FUNCTIONS_LIB = 
SOLVER          = 
NUMST           = 1
TID01EQ         = 0
NCSTATES        = 0
COMPUTER        = GLNX86
BUILDARGS       = ADD_MDL_NAME_TO_GLOBALS=1
MULTITASKING    = 0

LINUX_HOME = /usr/src/linux
RTAI_HOME  = /usr/src/rtai
COMEDI_HOME = /usr/src/comedi

#--------------------------- Tool Specifications -------------------------------

include $(MATLAB_ROOT)/rtw/c/tools/unixtools.mk

#------------------------------ Include Path -----------------------------------

MATLAB_INCLUDES = \
	-I$(MATLAB_ROOT)/simulink/include \
        -I$(MATLAB_ROOT)/extern/include \
	-I$(MATLAB_ROOT)/rtw/c/src \
	-I$(MATLAB_ROOT)/rtw/c/libsrc \
	-I$(MATLAB_ROOT)/rtw/c/rtai/devices
 
USER_INCLUDES = -I$(LINUX_HOME)/include -I$(RTAI_HOME)/include

# Additional file include paths
ADD_INCLUDES = \
       -I/usr/local/matlab/rtw/c/libsrc \


INCLUDES = -I. -I.. $(MATLAB_INCLUDES) $(ADD_INCLUDES) $(USER_INCLUDES) \
	$(INSTRUMENT_INCLUDES)

#-------------------------------- C Flags --------------------------------------

# Optimization Options
OPT_OPTS = -O2 

# General User Options
OPTS =
DEB_OPT = -DDBGPRT

# Compiler options, etc:
CC_OPTS = -Wall $(DEB_OPT) $(OPT_OPTS) $(OPTS) $(ANSI_OPTS) $(EXT_CC_OPTS) \
	$(LOG_OPTS)

CPP_REQ_DEFINES = -DMODEL=$(MODEL) -DRT -DNUMST=$(NUMST) \
		  -DTID01EQ=$(TID01EQ) -DNCSTATES=$(NCSTATES) -DUNIX \
		  -DMT=$(MULTITASKING)

CFLAGS = -DRT -DUSE_RTMODEL $(CC_OPTS) $(CPP_REQ_DEFINES) $(INCLUDES) $(NOVERSION) -ffast-math -c 

RT_MAIN_DEFINES = 
LDFLAGS = -static 

#-------------------------- Additional Libraries ------------------------------

SYSLIBS = $(EXT_LIB) -lpthread -lm

LIBS =

ifeq ($(OPT_OPTS),$(DEFAULT_OPT_OPTS))
LIBS += $(MATLAB_ROOT)/rtw/c/lib/glnx86/rtwlib_std.a
else
LIBS += rtwlib.a
endif
 
LIBS += $(S_FUNCTIONS_LIB) $(INSTRUMENT_LIBS)

#----------------------------- Source Files ------------------------------------

REQ_SRCS  = $(MODEL).c $(MODULES) rtmain.c rt_sim.c rt_nonfinite.c $(EXT_SRC)
USER_SRCS =
USER_OBJS = 
LOCAL_USER_OBJS = $(notdir $(USER_OBJS))

SRCS      = $(REQ_SRCS) $(S_FUNCTIONS) $(SOLVER)
OBJS      = $(SRCS:.c=.o) $(USER_OBJS)
LINK_OBJS = $(SRCS:.c=.o) $(LOCAL_USER_OBJS)

PROGRAM = ../$(MODEL)

#--------------------------------- Rules ---------------------------------------

$(PROGRAM) : $(OBJS) $(LIBS)
	gcc  $(LDFLAGS) -o $@ $(LINK_OBJS) $(LIBS) $(SYSLIBS)
	@echo "### Created executable: $(PROGRAM) ###"

%.o : %.c
	gcc -c  $(CFLAGS) $<

%.o : $(MATLAB_ROOT)/rtw/c/src/%.c
	gcc -c $(CFLAGS) $<

%.o : $(MATLAB_ROOT)/rtw/c/rtai/%.c
	gcc -c $(CFLAGS) $(RT_MAIN_DEFINES) $<

sfun_comedi%.o : $(MATLAB_ROOT)/rtw/c/rtai/devices/sfun_comedi%.c
	gcc -c $(CFLAGS) -I$(COMEDI_HOME)/include $<
	
sfun_mbdyn_com%.o : $(MATLAB_ROOT)/rtw/c/rtai/devices//sfun_mbdyn_com%.c
	gcc -c $(CFLAGS) $<
	
%.o : $(MATLAB_ROOT)/rtw/c/rtai/devices/%.c
	gcc -c $(CFLAGS) $<

%.o : $(MATLAB_ROOT)/rtw/c/libsrc/%.c
	gcc -c $(CFLAGS) $<

%.o : /usr/local/matlab/rtw/c/libsrc/%.c
	gcc -c $(CFLAGS) $<



%.o : $(MATLAB_ROOT)/simulink/src/%.c
	gcc -c $(CFLAGS) $<

#------------------------------- Libraries -------------------------------------



MODULES_rtwlib = \
    rt_atan2.o \
    rt_atan232.o \
    rt_backsubcc_dbl.o \
    rt_backsubcc_sgl.o \
    rt_backsubrc_dbl.o \
    rt_backsubrc_sgl.o \
    rt_backsubrr_dbl.o \
    rt_backsubrr_sgl.o \
    rt_enab.o \
    rt_forwardsubcc_dbl.o \
    rt_forwardsubcc_sgl.o \
    rt_forwardsubcr_dbl.o \
    rt_forwardsubcr_sgl.o \
    rt_forwardsubrc_dbl.o \
    rt_forwardsubrc_sgl.o \
    rt_forwardsubrr_dbl.o \
    rt_forwardsubrr_sgl.o \
    rt_hypot.o \
    rt_hypot32.o \
    rt_i32zcfcn.o \
    rt_intrp1flatcd.o \
    rt_intrp1flatcf.o \
    rt_intrp1flatd.o \
    rt_intrp1flatf.o \
    rt_intrp1lincd.o \
    rt_intrp1lincf.o \
    rt_intrp1linxd.o \
    rt_intrp1linxf.o \
    rt_intrp2flatcd.o \
    rt_intrp2flatcf.o \
    rt_intrp2flatd.o \
    rt_intrp2flatf.o \
    rt_intrp2lincd.o \
    rt_intrp2lincf.o \
    rt_intrp2linxd.o \
    rt_intrp2linxf.o \
    rt_intrpnflatcd.o \
    rt_intrpnflatcf.o \
    rt_intrpnflatd.o \
    rt_intrpnflatf.o \
    rt_intrpnlincd.o \
    rt_intrpnlincf.o \
    rt_intrpnlind.o \
    rt_intrpnlinf.o \
    rt_intrpnlinxd.o \
    rt_intrpnlinxf.o \
    rt_intrpnspld.o \
    rt_intrpnsplf.o \
    rt_logging.o \
    rt_look.o \
    rt_look1d.o \
    rt_look1d32.o \
    rt_look2d32_general.o \
    rt_look2d32_normal.o \
    rt_look2d_general.o \
    rt_look2d_normal.o \
    rt_look32.o \
    rt_lookflat1bincsd.o \
    rt_lookflat1bincsf.o \
    rt_lookflat1binczd.o \
    rt_lookflat1binczf.o \
    rt_lookflat1evncsd.o \
    rt_lookflat1evncsf.o \
    rt_lookflat1evnczd.o \
    rt_lookflat1evnczf.o \
    rt_lookflat1lincsd.o \
    rt_lookflat1lincsf.o \
    rt_lookflat1linczd.o \
    rt_lookflat1linczf.o \
    rt_lookflat2bincsd.o \
    rt_lookflat2bincsf.o \
    rt_lookflat2binczd.o \
    rt_lookflat2binczf.o \
    rt_lookflat2evncsd.o \
    rt_lookflat2evncsf.o \
    rt_lookflat2evnczd.o \
    rt_lookflat2evnczf.o \
    rt_lookflat2lincsd.o \
    rt_lookflat2lincsf.o \
    rt_lookflat2linczd.o \
    rt_lookflat2linczf.o \
    rt_lookflatnbincsd.o \
    rt_lookflatnbincsf.o \
    rt_lookflatnbinczd.o \
    rt_lookflatnbinczf.o \
    rt_lookflatnevncsd.o \
    rt_lookflatnevncsf.o \
    rt_lookflatnevnczd.o \
    rt_lookflatnevnczf.o \
    rt_lookflatnlincsd.o \
    rt_lookflatnlincsf.o \
    rt_lookflatnlinczd.o \
    rt_lookflatnlinczf.o \
    rt_looklin1bincsd.o \
    rt_looklin1bincsf.o \
    rt_looklin1binczd.o \
    rt_looklin1binczf.o \
    rt_looklin1binxsd.o \
    rt_looklin1binxsf.o \
    rt_looklin1binxzd.o \
    rt_looklin1binxzf.o \
    rt_looklin1evncsd.o \
    rt_looklin1evncsf.o \
    rt_looklin1evnczd.o \
    rt_looklin1evnczf.o \
    rt_looklin1evnxsd.o \
    rt_looklin1evnxsf.o \
    rt_looklin1evnxzd.o \
    rt_looklin1evnxzf.o \
    rt_looklin1lincsd.o \
    rt_looklin1lincsf.o \
    rt_looklin1linczd.o \
    rt_looklin1linczf.o \
    rt_looklin1linxsd.o \
    rt_looklin1linxsf.o \
    rt_looklin1linxzd.o \
    rt_looklin1linxzf.o \
    rt_looklin2bincsd.o \
    rt_looklin2bincsf.o \
    rt_looklin2binczd.o \
    rt_looklin2binczf.o \
    rt_looklin2binxsd.o \
    rt_looklin2binxsf.o \
    rt_looklin2binxzd.o \
    rt_looklin2binxzf.o \
    rt_looklin2evncsd.o \
    rt_looklin2evncsf.o \
    rt_looklin2evnczd.o \
    rt_looklin2evnczf.o \
    rt_looklin2evnxsd.o \
    rt_looklin2evnxsf.o \
    rt_looklin2evnxzd.o \
    rt_looklin2evnxzf.o \
    rt_looklin2lincsd.o \
    rt_looklin2lincsf.o \
    rt_looklin2linczd.o \
    rt_looklin2linczf.o \
    rt_looklin2linxsd.o \
    rt_looklin2linxsf.o \
    rt_looklin2linxzd.o \
    rt_looklin2linxzf.o \
    rt_looklinnbincsd.o \
    rt_looklinnbincsf.o \
    rt_looklinnbinczd.o \
    rt_looklinnbinczf.o \
    rt_looklinnbinxsd.o \
    rt_looklinnbinxsf.o \
    rt_looklinnbinxzd.o \
    rt_looklinnbinxzf.o \
    rt_looklinnevncsd.o \
    rt_looklinnevncsf.o \
    rt_looklinnevnczd.o \
    rt_looklinnevnczf.o \
    rt_looklinnevnxsd.o \
    rt_looklinnevnxsf.o \
    rt_looklinnevnxzd.o \
    rt_looklinnevnxzf.o \
    rt_looklinnlincsd.o \
    rt_looklinnlincsf.o \
    rt_looklinnlinczd.o \
    rt_looklinnlinczf.o \
    rt_looklinnlinxsd.o \
    rt_looklinnlinxsf.o \
    rt_looklinnlinxzd.o \
    rt_looklinnlinxzf.o \
    rt_looksplnbincsd.o \
    rt_looksplnbincsf.o \
    rt_looksplnbinczd.o \
    rt_looksplnbinczf.o \
    rt_looksplnbinssd.o \
    rt_looksplnbinssf.o \
    rt_looksplnbinszd.o \
    rt_looksplnbinszf.o \
    rt_looksplnbinxsd.o \
    rt_looksplnbinxsf.o \
    rt_looksplnbinxzd.o \
    rt_looksplnbinxzf.o \
    rt_looksplnevncsd.o \
    rt_looksplnevncsf.o \
    rt_looksplnevnczd.o \
    rt_looksplnevnczf.o \
    rt_looksplnevnssd.o \
    rt_looksplnevnssf.o \
    rt_looksplnevnszd.o \
    rt_looksplnevnszf.o \
    rt_looksplnevnxsd.o \
    rt_looksplnevnxsf.o \
    rt_looksplnevnxzd.o \
    rt_looksplnevnxzf.o \
    rt_looksplnlincsd.o \
    rt_looksplnlincsf.o \
    rt_looksplnlinczd.o \
    rt_looksplnlinczf.o \
    rt_looksplnlinssd.o \
    rt_looksplnlinssf.o \
    rt_looksplnlinszd.o \
    rt_looksplnlinszf.o \
    rt_looksplnlinxsd.o \
    rt_looksplnlinxsf.o \
    rt_looksplnlinxzd.o \
    rt_looksplnlinxzf.o \
    rt_lu_cplx.o \
    rt_lu_cplx_sgl.o \
    rt_lu_real.o \
    rt_lu_real_sgl.o \
    rt_matdivcc_dbl.o \
    rt_matdivcc_sgl.o \
    rt_matdivcr_dbl.o \
    rt_matdivcr_sgl.o \
    rt_matdivrc_dbl.o \
    rt_matdivrc_sgl.o \
    rt_matdivrr_dbl.o \
    rt_matdivrr_sgl.o \
    rt_matmultandinccc_dbl.o \
    rt_matmultandinccc_sgl.o \
    rt_matmultandinccr_dbl.o \
    rt_matmultandinccr_sgl.o \
    rt_matmultandincrc_dbl.o \
    rt_matmultandincrc_sgl.o \
    rt_matmultandincrr_dbl.o \
    rt_matmultandincrr_sgl.o \
    rt_matmultcc_dbl.o \
    rt_matmultcc_sgl.o \
    rt_matmultcr_dbl.o \
    rt_matmultcr_sgl.o \
    rt_matmultrc_dbl.o \
    rt_matmultrc_sgl.o \
    rt_matmultrr_dbl.o \
    rt_matmultrr_sgl.o \
    rt_matrx.o \
    rt_nrand.o \
    rt_plookbincd.o \
    rt_plookbincf.o \
    rt_plookbinkcd.o \
    rt_plookbinkcf.o \
    rt_plookbinxd.o \
    rt_plookbinxf.o \
    rt_plookevncd.o \
    rt_plookevncf.o \
    rt_plookevnkcd.o \
    rt_plookevnkcf.o \
    rt_plookevnxd.o \
    rt_plookevnxf.o \
    rt_plooklincd.o \
    rt_plooklincf.o \
    rt_plooklinkcd.o \
    rt_plooklinkcf.o \
    rt_plooklinxd.o \
    rt_plooklinxf.o \
    rt_printf.o \
    rt_sat_div_int16.o \
    rt_sat_div_int32.o \
    rt_sat_div_int8.o \
    rt_sat_div_uint16.o \
    rt_sat_div_uint32.o \
    rt_sat_div_uint8.o \
    rt_sat_prod_int16.o \
    rt_sat_prod_int32.o \
    rt_sat_prod_int8.o \
    rt_sat_prod_uint16.o \
    rt_sat_prod_uint32.o \
    rt_sat_prod_uint8.o \
    rt_tdelay.o \
    rt_urand.o \
    rt_zcfcn.o \


rtwlib.a : $(MAKEFILE) rtw_proj.tmw $(MODULES_rtwlib)
	@echo "### Creating $@ "
	$(AR) r $@ $(MODULES_rtwlib)



#----------------------------- Dependencies ------------------------------------

$(OBJS) : $(MAKEFILE) rtw_proj.tmw

#--------- Miscellaneous rules to purge and clean ------------------------------

purge : clean
	@echo "### Deleting the generated source code for $(MODEL)"
	@\rm -f $(MODEL).c $(MODEL).h $(MODEL).prm $(MODEL).reg $(MODEL).rtw \
	        $(MODULES) rtw_proj.tmw $(MAKEFILE)

clean :
	@echo "### Deleting the objects and $(PROGRAM)"
	@\rm -f $(MODULES_rtwlib) $(MODULES_dsp_rt)
	@\rm -f $(LINK_OBJS) $(PROGRAM)

# EOF: pendulum_control.mk
