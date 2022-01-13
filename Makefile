all emcc clean:
	$(MAKE) -C libsrc $@

debug:
	EMCC_DEBUG=1 $(MAKE) -C libsrc $@

FINDCMD = find trimesh2 -name 'OBJ*' -prune -o -name '.git*' -prune -o -type f -print

tar:
	cd .. && tar zcvf trimesh2.tar.gz `$(FINDCMD) | sort`

zip:
	cd .. && $(FINDCMD) | sort | zip -9 trimesh2 -@

.PHONY : all clean debug tar zip emcc
