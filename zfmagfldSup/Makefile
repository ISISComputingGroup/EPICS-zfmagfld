TOP=..
include $(TOP)/configure/CONFIG
#=======================================

# Install .dbd and .db files
DB += zfmagfld.db
DB += zfmagfld_axes.db
DB += zfmagfld_sensor_matrix.db
DB += zfmagfld_errors.db
DB += zfmagfld_extra_axis.db

LIBRARY_IOC = ZFMAGFLD
DBD += zfmagfld.dbd

ZFMAGFLD_DBD += asyn.dbd
ZFMAGFLD_DBD += aSubFunctions.dbd

ZFMAGFLD_SRCS += zfmagfld.c zfmagfld_impl.cpp

ZFMAGFLD_LIBS += $(EPICS_BASE_IOC_LIBS)
ZFMAGFLD_LIBS += asyn utilities gsl
ZFMAGFLD_LIBS += asubFunctions

#=======================================
include $(TOP)/configure/RULES
