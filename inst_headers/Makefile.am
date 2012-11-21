xtincludedir = $(includedir)/xt

include_HEADERS = yaxt.h

BUILT_SOURCES = yaxt.h

SUBDIRS = f90

yaxt.h: Makefile
	( echo '/* license */' ; echo '#ifndef YAXT_H' ; echo '#define YAXT_H' ; \
	for header in $(xtinclude_HEADERS) ; \
	  do echo "#include <$$header>" ; done ; \
	echo '#endif' \
	) >yaxt.h.tmp
	if test ! -f $@ || ! cmp $@ $@.tmp >/dev/null ; then \
	mv $@.tmp $@ ; else rm $@.tmp ; fi

all-local:
	$(MKDIR_P) xt
	for header in $(xtinclude_HEADERS) ; \
	  do \
	  headerbasedir=`test -f '../src/$$header' || echo '$(srcdir)/'` ; \
          if test "$$headerbasedir../src/$$header" -ot "$$header" \
	    -a ../include/config.h -ot "$$header" ; then : ; else \
	  $(PERL) -I$(top_builddir)/scripts \
	  $(top_builddir)/scripts/header2installheader.pl \
	  --config-header=../include/config.h \
	  --srcdir=$$headerbasedir../src \
	  --dstdir=. \
	  $$headerbasedir../src/$$header ; \
	  fi ; \
	  done

clean-local:
	rm -rf xt yaxt.h

install-data-local:
	@echo " $(MKDIR_P) '$(DESTDIR)$(xtincludedir)'"; \
          $(MKDIR_P) "$(DESTDIR)$(xtincludedir)" || exit 1; \
	for p in $(xtinclude_HEADERS); do \
          echo " $(INSTALL_HEADER) $$p '$(DESTDIR)$(xtincludedir)'"; \
          $(INSTALL_HEADER) $$p "$(DESTDIR)$(xtincludedir)" || exit $$?; \
        done

uninstall-local:
	rm -rf "$(DESTDIR)$(xtincludedir)"

dist-hook:
	-rm $(distdir)/xt/*.h
	-rmdir $(distdir)/xt

@AMDEP_TRUE@@am__include@ @am__quote@../src/xtinclude_Makefile.inc@am__quote@
