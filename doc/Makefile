
doc: Doxyfile
	cd .. && doxygen doc/Doxyfile

Doxyfile: Doxyfile.in ../configure.ac
	./mk_doxyfile.pl

clean:
	rm -rf html/

