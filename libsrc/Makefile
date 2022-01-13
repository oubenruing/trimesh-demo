DESTDIR = ../dist.$(UNAME)
INCLUDES = -I../include

include ../Makerules

CCFILES =	TriMesh.cc \
		TriMesh_bounding.cc \
		TriMesh_connectivity.cc \
		TriMesh_curvature.cc \
		TriMesh_io.cc \
		TriMesh_grid.cc \
		TriMesh_normals.cc \
		TriMesh_pointareas.cc \
		TriMesh_stats.cc \
		TriMesh_tstrips.cc \
		GLCamera.cc \
		GLManager.cc \
		ICP.cc \
		KDtree.cc \
		conn_comps.cc \
		diffuse.cc \
		edgeflip.cc \
		faceflip.cc \
		filter.cc \
		make.cc \
		merge.cc \
		overlap.cc \
		remove.cc \
		reorder_verts.cc \
		subdiv.cc \
		umbrella.cc \



OFILES = $(addprefix $(OBJDIR)/,$(CCFILES:.cc=.o))
OUTPUT = $(DESTDIR)/libtrimesh.html
default: $(OUTPUT)

$(OUTPUT): $(OFILES)
	$(EMCC)

clean:
	-rm -f $(OBJDIR)/*.o $(OBJDIR)/*.d
	-rmdir $(OBJDIR)

spotless: clean
	-rm -f $(LIB)
	-rmdir $(DESTDIR)

