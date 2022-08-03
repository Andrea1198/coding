# Copyright (C) 2022-20 Sprel_sium
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License. See the file `License' in the root directory
# of the present distribution.

include make.inc

# execute a target irrespective of the presence of a file or directory
# with the same name
.PHONY: install

default :
        @echo 'to install Quantum ESPRESSO, type at the shell prompt:'
        @echo '  ./configure [--prefix=]'
        @echo '  make [-j] target'
        @echo ' '
        @echo 'where target identifies one or multiple CORE PACKAGES:'
        @echo '  pw           basic code for scf, structure optimization, MD'
        @echo '  ph           phonon code, Gamma-only and third-order derivatives'
        @echo '  hp           calculation of the Hubbard parameters from DFPT'
        @echo '  pwcond       ballistic conductance'
        @echo '  neb          code for Nudged Elastic Band method'
        @echo '  pp           postprocessing programs'
        @echo '  pwall        same as "make pw ph pp pwcond neb"'
        @echo '  cp           CP code: Car-Parrinello molecular dynamics'
        @echo '  all_currents QEHeat code: energy flux and charge current'
        @echo '  tddfpt       time dependent dft code'
        @echo '  gwl          GW with Lanczos chains'
        @echo '  ld1          utilities for pseudopotential generation'
        @echo '  xspectra     X-ray core-hole spectroscopy calculations'
        @echo '  couple       Library interface for coupling to external codes'
        @echo '  epw          Electron-Phonon Coupling with Wannier functions'
        @echo '  gui          Graphical User Interface'
        @echo '  all          same as "make pwall cp ld1 tddfpt xspectra hp"'
        @echo ' '
        @echo 'where target identifies one or multiple THIRD-PARTIES PACKAGES:'
        @echo '  gipaw        NMR and EPR spectra'
        @echo '  w90          Maximally localised Wannier Functions'
        @echo '  want         Quantum Transport with Wannier functions'
        @echo '  yambo        electronic excitations with plane waves'
        @echo '  d3q          general third-order code and thermal transport codes'
        @echo ' '
        @echo 'where target is one of the following suite operation:'
        @echo '  doc          build documentation'
        @echo '  links        create links to all executables in bin/'
        @echo '  install      copy all executables to PREFIX/bin/'
        @echo '               (works with "configure --prefix=PREFIX)"'
        @echo '  tar          create a tarball of the source tree'
        @echo '  depend       generate dependencies (make.depend files)'
        @if test -d GUI/; then \
                echo '  tar-gui      create a standalone PWgui tarball from the GUI sources'; \
                echo '  tar-qe-modes create a tarball for QE-modes (Emacs major modes for Quantum ESPRESSO)'; fi
        @echo '  clean        remove executables and objects'
        @echo '  veryclean    remove files produced by "configure" as well'
        @echo '  distclean    revert distribution to the original status'
