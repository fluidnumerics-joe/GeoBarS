
SHELL=/bin/bash
GeoBarSDIR=$(shell pwd | sed 's/\GeoBarS.*/GeoBarS/')
INITDIR=$(shell pwd)

.PHONY: geobars

geobars :
	@echo 'Found installation directory : '${GeoBarSDIR}
	make MODDIR=${INITDIR} -C ${GeoBarSDIR}/build/ geobars
	mv ${GeoBarSDIR}/build/geobars ./

clean :
	make -C ${GeoBarSDIR}/build/ clean

