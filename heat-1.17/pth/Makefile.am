noinst_PROGRAMS = heat barrier pipeline histo

heat_SOURCES = $(COMMON)/main.c compute.c
barrier_SOURCES = barrier.c
pipeline_SOURCES = pipeline.c
histo_SOURCES = histo.c

## Uncomment the following:
# noinst_PROGRAMS += barrier_test
# barrier_test_SOURCES = ...

#
# Maintenance rules
#
CONFIGURE_DEPENDENCIES = $(srcdir)/../VERSION
ACLOCAL_AMFLAGS = --force -I .
include $(srcdir)/../src/libcommon.mk


