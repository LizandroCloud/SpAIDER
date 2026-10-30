# Copyright 1994-2008 The MathWorks, Inc.
#
# File    : raccel_vc.tmf   $Revision: 1.1.6.5.4.1 $
#
# Abstract:
#       Real-Time Workshop template makefile for building a PC-based "rapid
#       acceleration" executable from the generated C code for a Simulink model
#       using generated C code and the
#               Microsoft Visual C/C++ compiler versions: 6.0, 7.1, 8.0
#
#       Note that this template is automatically customized by the Real-Time
#       Workshop build procedure to create "<model>.mk"
#
#       The following defines can be used to modify the behavior of the
#       build:
#
#         OPT_OPTS       - Optimization option. See DEFAULT_OPT_OPTS in
#                          vctools.mak for default.
#         OPTS           - User specific options.
#         CPP_OPTS       - C++ compiler options.
#         USER_SRCS      - Additional user sources, such as files needed by
#                          S-functions.
#         USER_INCLUDES  - Additional include paths
#                          (i.e. USER_INCLUDES="-Iwhere-ever -Iwhere-ever2")
#
#       To enable debugging:
#         set OPT_OPTS=-Zi (may vary with compiler version, see compiler doc) 
#         modify tmf LDFLAGS to include /DEBUG
#
#       This template makefile is designed to be used with a system target
#       file that contains 'rtwgensettings.BuildDirSuffix' see raccel.tlc


#------------------------ Macros read by make_rtw -----------------------------
#
# The following macros are read by the Real-Time Workshop build procedure:
#
#  MAKECMD         - This is the command used to invoke the make utility
#  HOST            - What platform this template makefile is targeted for
#                    (i.e. PC or UNIX)
#  BUILD           - Invoke make from the Real-Time Workshop build procedure
#                    (yes/no)?
#  SYS_TARGET_FILE - Name of system target file.

MAKECMD         = nmake
HOST            = PC
BUILD           = yes
SYS_TARGET_FILE = raccel.tlc
BUILD_SUCCESS	= *** Created
COMPILER_TOOL_CHAIN = vc

#---------------------- Tokens expanded by make_rtw ---------------------------
#
# The following tokens, when wrapped with "|>" and "<|" are expanded by the
# Real-Time Workshop build procedure.
#
#  MODEL_NAME          - Name of the Simulink block diagram
#  MODEL_MODULES       - Any additional generated source modules
#  MAKEFILE_NAME       - Name of makefile created from template makefile <model>.mk
#  MATLAB_ROOT         - Path to were MATLAB is installed.
#  MATLAB_BIN          - Path to MATLAB executable.
#  S_FUNCTIONS         - List of S-functions.
#  S_FUNCTIONS_LIB     - List of S-functions libraries to link.
#  BUILDARGS           - Options passed in at the command line.
#  EXT_MODE            - yes (1) or no (0): Build for external mode
#  TMW_EXTMODE_TESTING - yes (1) or no (0): Build ext_test.c for external mode
#                        testing.
#  EXTMODE_TRANSPORT   - Name of transport mechanism (e.g. tcpip, serial) for extmode
#  EXTMODE_STATIC      - yes (1) or no (0): Use static instead of dynamic mem alloc.
#  EXTMODE_STATIC_SIZE - Size of static memory allocation buffer.
#  MEX_OPT_FILE        - mex option file name
#  VISUAL_VER          - Visual C++ version information


MODEL           = emso2
MODULES         = emso2_data.c rt_logging.c rt_matrx.c rt_nonfinite.c 
MAKEFILE        = emso2.mk
MATLAB_ROOT     = E:\Program Files\MATLAB2008\R2008a
ALT_MATLAB_ROOT = E:\PROGRA~1\MATLAB2008\R2008a
MATLAB_BIN      = E:\Program Files\MATLAB2008\R2008a\bin
ALT_MATLAB_BIN  = E:\PROGRA~1\MATLAB2008\R2008a\bin
S_FUNCTIONS     = emso_sf_c.c
S_FUNCTIONS_LIB = 
BUILDARGS       =  GENERATE_REPORT=0 RSIM_SOLVER_SELECTION=2
RSIM_PARAMETER_LOADING = 1
MEX_OPT_FILE    = -f E:\DOCUME~1\ADMINI~1\APPLIC~1\MATHWO~1\MATLAB\R2008a\mexopts.bat
VISUAL_VER      = 6.0

EXT_MODE            = 1
TMW_EXTMODE_TESTING = 0
EXTMODE_TRANSPORT   = 0
EXTMODE_STATIC      = 0
EXTMODE_STATIC_SIZE = 1000000

SOLVER              = ode23s.c
SOLVER_TYPE         = VariableStep
NUMST               = 2
TID01EQ             = 0
NCSTATES            = 40
MULTITASKING        = 0
PCMATLABROOT        = E:\\Program Files\\MATLAB2008\\R2008a

MODELREFS            = 
SHARED_SRC           = 
SHARED_SRC_DIR       = 
SHARED_BIN_DIR       = 
SHARED_LIB           = 
GEN_SAMPLE_MAIN      = 0

OPTIMIZATION_FLAGS   = /Od /Oy- /DNDEBUG
ADDITIONAL_LDFLAGS   = 

#--------------------------- Model and reference models -----------------------
MODELLIB                  = emso2lib.lib
MODELREF_LINK_LIBS        = 
MODELREF_LINK_RSPFILE     = emso2_ref.rsp
MODELREF_INC_PATH         = 
RELATIVE_PATH_TO_ANCHOR   = ..
MODELREF_TARGET_TYPE      = NONE

GLOBAL_TIMING_ENGINE       = 0

!if "$(MATLAB_ROOT)" != "$(ALT_MATLAB_ROOT)"
MATLAB_ROOT = $(ALT_MATLAB_ROOT)
!endif
!if "$(MATLAB_BIN)" != "$(ALT_MATLAB_BIN)"
MATLAB_BIN = $(ALT_MATLAB_BIN)
!endif

#---------------------------Solver---------------------------------------------

RSIM_WITH_SL_SOLVER = 1

#--------------------------- Tool Specifications ------------------------------
!include $(MATLAB_ROOT)\rtw\c\tools\vctools.mak

PERL = $(MATLAB_ROOT)\sys\perl\win32\bin\perl

#------------------------------ Include Path ----------------------------------

MATLAB_INCLUDES =                    $(MATLAB_ROOT)\simulink\include
MATLAB_INCLUDES = $(MATLAB_INCLUDES);$(MATLAB_ROOT)\extern\include
MATLAB_INCLUDES = $(MATLAB_INCLUDES);$(MATLAB_ROOT)\rtw\c\src
MATLAB_INCLUDES = $(MATLAB_INCLUDES);$(MATLAB_ROOT)\rtw\c\raccel
MATLAB_INCLUDES = $(MATLAB_INCLUDES);$(MATLAB_ROOT)\rtw\c\src\ext_mode\common

# Additional includes

MATLAB_INCLUDES = $(MATLAB_INCLUDES);C:\Dropbox\Spaider\problems\emso\AceticAnhydride_Process\emso2_raccel_rtw

MATLAB_INCLUDES = $(MATLAB_INCLUDES);C:\Dropbox\Spaider\problems\emso\AceticAnhydride_Process

MATLAB_INCLUDES = $(MATLAB_INCLUDES);$(MATLAB_ROOT)\rtw\c\libsrc


INCLUDE = .;$(RELATIVE_PATH_TO_ANCHOR);$(MATLAB_INCLUDES);$(INCLUDE)

!if "$(SHARED_SRC_DIR)" != ""
INCLUDE = $(INCLUDE);$(SHARED_SRC_DIR)
!endif

#------------------------ External mode ---------------------------------------
# Uncomment -DVERBOSE to have information printed to stdout
# To add a new transport layer, see the comments in
#   <matlabroot>/toolbox/simulink/simulink/extmode_transports.m

EXT_CC_OPTS = -DEXT_MODE # -DVERBOSE
!if $(EXTMODE_TRANSPORT) == 0 #tcpip
EXT_SRC = ext_svr.c updown.c ext_work.c rtiostream_interface.c rtiostream_tcpip.c
EXT_LIB = wsock32.lib
!endif
!if $(EXTMODE_TRANSPORT) == 1 #serial_win32
EXT_SRC = ext_svr.c updown.c ext_work.c ext_svr_serial_transport.c
EXT_SRC = $(EXT_SRC) ext_serial_pkt.c ext_serial_win32_port.c
EXT_LIB =
!endif
!if $(TMW_EXTMODE_TESTING) == 1
EXT_SRC     = $(EXT_SRC) ext_test.c
EXT_CC_OPTS = $(EXT_CC_OPTS) -DTMW_EXTMODE_TESTING
!endif
!if $(EXTMODE_STATIC) == 1
EXT_SRC     = $(EXT_SRC) mem_mgr.c
EXT_CC_OPTS = $(EXT_CC_OPTS) -DEXTMODE_STATIC -DEXTMODE_STATIC_SIZE=$(EXTMODE_STATIC_SIZE)
!endif

#----------------------------- Parameter Tuning -------------------------------
PARAM_CC_OPTS = -DRSIM_PARAMETER_LOADING

#----------------- Compiler and Linker Options --------------------------------

# Optimization Options
#   Set OPT_OPTS=-Zi for debugging (may depend on compiler version)
OPT_OPTS = $(DEFAULT_OPT_OPTS)

# General User Options
OPTS =

!if "$(OPTIMIZATION_FLAGS)" != ""
CC_OPTS = $(OPTS) $(EXT_CC_OPTS) $(OPTIMIZATION_FLAGS)
!else
CC_OPTS = $(OPT_OPTS) $(OPTS) $(EXT_CC_OPTS)
!endif


CPP_REQ_DEFINES = -DMODEL=$(MODEL) -DHAVESTDIO

CPP_REQ_DEFINES = $(CPP_REQ_DEFINES) \
                  -DNRT -DRSIM_WITH_SL_SOLVER -DMATLAB_ROOT="$(PCMATLABROOT)"
!if "$(MULTITASKING)" == "1"
CPP_REQ_DEFINES = $(CPP_REQ_DEFINES) \
                  -DRSIM_WITH_SOLVER_MULTITASKING \
                  -DTID01EQ=$(TID01EQ) \
                  -DNUMST=$(NUMST)
!endif

# Add flags to generate the pdb file if -DGEN_PDB=1 is specified
GEN_PDB = 0
!if $(GEN_PDB) == 1
CPP_REQ_DEFINES = $(CPP_REQ_DEFINES) /Zi /Fd"$(MODEL).pdb"
!endif

# Uncomment this line to move warning level to W4
# cflags = $(cflags:W3=W4)
CFLAGS   = $(MODELREF_INC_PATH) $(cflags) $(cvarsmt) $(PARAM_CC_OPTS) \
	   $(CFLAGS_VERSPECIFIC) $(CC_OPTS) $(CPP_REQ_DEFINES) $(USER_INCLUDES)
CPPFLAGS = $(MODELREF_INC_PATH) $(cflags) $(cvarsmt) $(PARAM_CC_OPTS) \
	   $(CPPFLAGS_VERSPECIFIC) $(CPP_OPTS) $(CC_OPTS) $(CPP_REQ_DEFINES) \
	   $(USER_INCLUDES)

!if $(GEN_PDB) == 1
LDFLAGS  = /DEBUG /PDB:"$(MODEL).pdb"
!else
LDFLAGS  = $(ldebug)
!endif
LDFLAGS  = $(LDFLAGS) $(conflags) $(EXT_LIB) $(conlibs) $(ADDITIONAL_LDFLAGS)

#----------------------------- Source Files -----------------------------------

#Standalone executable
!if "$(MODELREF_TARGET_TYPE)" == "NONE"
PRODUCT            = $(MODEL).exe	
BUILD_PRODUCT_TYPE = executable
REQ_SRCS           = $(MODEL).c $(MODULES) \
                     raccel_main.c raccel_sup.c raccel_mat.c \
                     $(EXT_SRC)
#Model Reference Target
!else
PRODUCT            = $(MODELLIB)
BUILD_PRODUCT_TYPE = library
REQ_SRCS           = $(MODULES) $(EXT_SRC)
!endif

USER_SRCS =

SRCS = $(REQ_SRCS) $(USER_SRCS) $(S_FUNCTIONS)
OBJS_CPP_UPPER = $(SRCS:.CPP=.obj)
OBJS_CPP_LOWER = $(OBJS_CPP_UPPER:.cpp=.obj)
OBJS_C_UPPER = $(OBJS_CPP_LOWER:.C=.obj)
OBJS = $(OBJS_C_UPPER:.c=.obj)
SHARED_OBJS = $(SHARED_SRC:.c=.obj)

#--------------------------- Required Libraries -------------------------------

MAT_LIBS = $(MATLAB_ROOT)\extern\lib\win32\microsoft\libut.lib \
           $(MATLAB_ROOT)\extern\lib\win32\microsoft\libmx.lib \
           $(MATLAB_ROOT)\extern\lib\win32\microsoft\libmat.lib \
	   $(MATLAB_ROOT)\extern\lib\win32\microsoft\libmwmathutil.lib \
           $(MATLAB_ROOT)\rtw\c\lib\win32\rsim_engine.lib

# ------------------------- Additional Libraries ------------------------------

LIBS =


LIBS = $(LIBS)

# ---------------------------- Linker Script ----------------------------------

CMD_FILE = $(MODEL).lnk
GEN_LNK_SCRIPT = $(MATLAB_ROOT)\rtw\c\tools\mkvc_lnk.pl

#--------------------------------- Rules --------------------------------------
all: set_environment_variables $(PRODUCT)

!if "$(MODELREF_TARGET_TYPE)" == "NONE"
#--- Stand-alone model ---
$(PRODUCT) : $(OBJS) $(SHARED_LIB) $(LIBS) $(MODELREF_LINK_LIBS)
	$(PERL) $(GEN_LNK_SCRIPT) $(CMD_FILE) $(OBJS)
	$(LD) $(LDFLAGS) $(S_FUNCTIONS_LIB) $(SHARED_LIB) $(LIBS) $(MAT_LIBS) \
		@$(CMD_FILE) @$(MODELREF_LINK_RSPFILE) -out:$@
	@del $(CMD_FILE)
!else
#--- Model reference RTW Target ---
$(PRODUCT) : $(OBJS) $(SHARED_LIB)
	$(PERL) $(GEN_LNK_SCRIPT) $(CMD_FILE) $(OBJS)
	$(LD) -lib /OUT:$@ @$(CMD_FILE) $(S_FUNCTIONS_LIB)
!endif
	@cmd /C "echo $(BUILD_SUCCESS) $(BUILD_PRODUCT_TYPE): $@"

!if $(GEN_SAMPLE_MAIN) == 0
{$(MATLAB_ROOT)\rtw\c\raccel}.c.obj :
	$(CC) $(CFLAGS) $<
!endif

{$(MATLAB_ROOT)\rtw\c\src}.c.obj :
	$(CC) $(CFLAGS) $<

{$(MATLAB_ROOT)\rtw\c\src\ext_mode\common}.c.obj :
	$(CC) $(CFLAGS) $<

{$(MATLAB_ROOT)\rtw\c\src\ext_mode\tcpip}.c.obj :
	$(CC) $(CFLAGS) $<

{$(MATLAB_ROOT)\rtw\c\src\ext_mode\serial}.c.obj :
	$(CC) $(CFLAGS) $<

{$(MATLAB_ROOT)\rtw\c\src\ext_mode\custom}.c.obj :
	$(CC) $(CFLAGS) $<

# Additional sources

{$(MATLAB_ROOT)\rtw\c\src}.c.obj :
	$(CC) $(CFLAGS) $<



{$(MATLAB_ROOT)\rtw\c\src}.cpp.obj :
	$(CC) $(CPPFLAGS) $<



# Look in simulink/src helper files

{$(MATLAB_ROOT)\simulink\src}.c.obj :
	$(CC) $(CFLAGS) $<

# Put these rule last, otherwise nmake will check toolboxes first

{$(RELATIVE_PATH_TO_ANCHOR)}.cpp.obj :
	$(CC) $(CPPFLAGS) $<

.cpp.obj :
	$(CC) $(CPPFLAGS) $<

{$(RELATIVE_PATH_TO_ANCHOR)}.c.obj :
	$(CC) $(CFLAGS) $<

.c.obj :
	$(CC) $(CFLAGS) $<

!if "$(SHARED_LIB)" != ""
$(SHARED_LIB) : $(SHARED_SRC)
	@$(CC) $(CFLAGS) -Fo$(SHARED_BIN_DIR)\ @<<
$?
<<
	@$(LIBCMD) /nologo /out:$@ $(SHARED_OBJS)
!endif

set_environment_variables:
	@set INCLUDE=$(INCLUDE)
	@set LIB=$(LIB)


# Libraries:





#----------------------------- Dependencies -----------------------------------

$(OBJS) : $(MAKEFILE) rtw_proj.tmw
