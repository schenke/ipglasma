all	: ipglasma

ipglasma:
	$(MAKE) -C src
	cp -f src/ipglasma ./

clean:
	$(MAKE) -C src clean
	rm -f ipglasma

distclean:
	$(MAKE) -C src distclean
	rm -f ipglasma
