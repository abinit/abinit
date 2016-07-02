#
# Makefile for the ABINIT Configurator
#

# Project name
PROJECT = abinit_config

# PyQT resource compiler
PYRCC       = pyrcc4
PYRCC_FLAGS =
PYRCC_NAME  = abinit_config_init
RESOURCES   =

# PyQT UI form compiler
PYUIC       = pyuic4
PYUIC_FLAGS = -x -i 2


# Targets
all_targets all: ui_$(PROJECT).py

ui_$(PROJECT).py: $(PROJECT).ui
	$(PYUIC) $(PYUIC_FLAGS) -o $@ $<

pyqt_resources:
	$(PYRCC) $(PYRCC_FLAGS) -name $(PYRCC_NAME) -o $(PROJECT)_rc.py $(RESOURCES)

clean:
	rm -f *.pyc *.pyo
	rm -f ui_$(PROJECT).py $(PROJECT)_rc.py
