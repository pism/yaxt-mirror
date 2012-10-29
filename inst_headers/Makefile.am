xtincludedir = $(includedir)/xt

include_HEADERS = yaxt.h

yaxt.h: Makefile
	for header in $(xtinclude_HEADERS) ; \
	  do echo "#include <$$header>" ; done \
	  >yaxt.h

all-local:
	$(MKDIR_P) xt
	for header in $(xtinclude_HEADERS) ; \
	  do ln -sf \
	  ../`test -f '../src/$$header' || echo '$(srcdir)/'`../src/$$header \
	  $$header ; done

install-data-local:
	@echo " $(MKDIR_P) '$(DESTDIR)$(xtincludedir)'"; \
          $(MKDIR_P) "$(DESTDIR)$(xtincludedir)" || exit 1; \
	for p in $(xtinclude_HEADERS); do \
          echo " $(INSTALL_HEADER) $$p '$(DESTDIR)$(xtincludedir)'"; \
          $(INSTALL_HEADER) $$p "$(DESTDIR)$(xtincludedir)" || exit $$?; \
        done

dist-hook:
	rm xt/*.h
	rmdir $(distdir)/xt

@AMDEP_TRUE@@am__include@ @am__quote@../src/xtinclude_Makefile.inc@am__quote@
