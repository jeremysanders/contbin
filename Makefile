# Makefile for contbin, Jeremy Sanders
######################################
# 

# add -lsocket to below for Solaris
# add -Ldirname to add directory to link path to look for cfitsio
linkflags=-lcfitsio -Lparammm -lparammm

# where to install (not very well tested)
bindir=/usr/local/bin

# sensible compiler flags
export CXXFLAGS=-O2 -g -Wall -std=c++11
export CXX=g++

##############################################################################
# you probably don't need to change anything below here

default: all

programs=contbin accumulate_smooth accumulate_smooth_expmap make_region_files \
	paint_output_images dumpdata exposure_smooth accumulate_smooth_expcorr \
	adaptive_gaussian_smooth
all: $(programs)

clean:
	-rm -f *.o parammm/*.o parammm/*.a $(programs)

install:
	install $(programs) $(bindir)

exposure_smooth.o: exposure_smooth.cc
adaptive_gaussian_smooth.o: adaptive_gaussian_smooth.cc
dumpdata.o: dumpdata.cc
binner.o: point.hh binner.cc binner.hh misc.hh bin.hh \
	scrubber.hh terminal.hh
contbin.o: binner.hh contbin.cc misc.hh
flux_estimator.o: flux_estimator.cc misc.hh flux_estimator.hh
bin.o: bin.hh bin.cc
scrubber.o: scrubber.cc scrubber.hh bin.hh
terminal.o: terminal.hh terminal.cc

adaptive_gaussian_smooth_objs=adaptive_gaussian_smooth.o fitsio_simple.o memimage.o

adaptive_gaussian_smooth: $(adaptive_gaussian_smooth_objs)  parammm/libparammm.a
	$(CXX) -o adaptive_gaussian_smooth $(adaptive_gaussian_smooth_objs) $(linkflags)

exposure_smooth_objs=exposure_smooth.o fitsio_simple.o memimage.o

exposure_smooth: $(exposure_smooth_objs)  parammm/libparammm.a
	$(CXX) -o exposure_smooth $(exposure_smooth_objs) $(linkflags)

dumpdata: dumpdata.o fitsio_simple.o
	$(CXX) -o dumpdata fitsio_simple.o dumpdata.o $(linkflags)

contbin_objs=contbin.o binner.o flux_estimator.o bin.o scrubber.o \
	terminal.o fitsio_simple.o memimage.o

contbin: $(contbin_objs) parammm/libparammm.a
	$(CXX) -o contbin $(contbin_objs) $(linkflags)

acc_smooth_objs=accumulate_smooth.o flux_estimator.o \
	fitsio_simple.o memimage.o

accumulate_smooth: $(acc_smooth_objs) parammm/libparammm.a
	$(CXX) -o accumulate_smooth $(acc_smooth_objs) $(linkflags)

acc_smooth_expmap_objs=accumulate_smooth_expmap.o \
	fitsio_simple.o memimage.o

accumulate_smooth_expmap: $(acc_smooth_expmap_objs) parammm/libparammm.a
	$(CXX) -o accumulate_smooth_expmap $(acc_smooth_expmap_objs) \
		$(linkflags)

acc_smooth_expcorr_objs=accumulate_smooth_expcorr.o \
	fitsio_simple.o memimage.o

accumulate_smooth_expcorr: $(acc_smooth_expcorr_objs) parammm/libparammm.a
	$(CXX) -o accumulate_smooth_expcorr $(acc_smooth_expcorr_objs) \
		$(linkflags)

make_region_files_objs=make_region_files.o \
        fitsio_simple.o memimage.o

make_region_files: $(make_region_files_objs) parammm/libparammm.a
	$(CXX) -o make_region_files $(make_region_files_objs) \
		$(linkflags)

paint_output_images_objs=paint_output_images.o \
	fitsio_simple.o memimage.o format_string.o

paint_output_images: $(paint_output_images_objs) parammm/libparammm.a
	$(CXX) -o paint_output_images $(paint_output_images_objs) \
		$(linkflags)

parammm/libparammm.a:
	$(MAKE) -C parammm libparammm.a
