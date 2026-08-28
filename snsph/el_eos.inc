c © 2026. Triad National Security, LLC. All rights reserved.
c This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos National Laboratory (LANL), which is operated by Triad National Security, LLC for the U.S. Department of Energy/National Nuclear Security Administration. All rights in the program are reserved by Triad National Security, LLC, and the U.S. Department of Energy/National Nuclear Security Administration. The Government is granted for itself and others acting on its behalf a nonexclusive, paid-up, irrevocable worldwide license in this material to reproduce, prepare. derivative works, distribute copies to the public, perform publicly and display publicly, and to permit others to do so.

C23456789012345678901234567890123456789012345678901234567890123456789012
C***********************************************************************
C
C    MODULE:       EL_EOS.INC
C    TYPE:         INCLUDE FILE
C    AUTHOR:       F. DOUGLAS SWESTY
C    DATE:         2/12/91
C
C
C    CALL LINE:    INCLUDE 'EL_EOS.INC'
C
C***********************************************************************
C23456789012345678901234567890123456789012345678901234567890123456789012
C
C
C                              Electron density and concentration
      DOUBLE PRECISION NSUBE, QSUBE
C
C                              Temporary coefficients
      DOUBLE PRECISION ACOEF, BCOEF, CCOEF, DBDT
C
C                              Electron chemical potential
      DOUBLE PRECISION MUSUBE
C
C                              Electron pressure, internal energy,
C                              Helmholtz free energy, and entropy
      DOUBLE PRECISION EPRESS, EU, FSUBE, ES
C
C                              Photon pressure, internal energy,
C                              Helmholtz free energy, and entropy
      DOUBLE PRECISION PPRESS, PS, PU, PF
C
C
      DOUBLE PRECISION DEMUDT, DEMUDN, DEMUDY
      DOUBLE PRECISION DEPDT, DEPDN, DEPDY
      DOUBLE PRECISION DESDT, DESDN, DESDY
      DOUBLE PRECISION DEUDT, DEUDN, DEUDY
C
C
      DOUBLE PRECISION DPPDN, DPPDT, DPPDY
      DOUBLE PRECISION DPSDN, DPSDT, DPSDY
      DOUBLE PRECISION DPUDN, DPUDT, DPUDY
C
C
C                              This common block contains all the
C                              electron and phton variables needed
C                              by the outside world
      COMMON /ELECTR/ NSUBE, MUSUBE, 
     1  EPRESS, EU, FSUBE, ES, PPRESS, PS, PU, PF,
     2  DEMUDT, DEMUDN, DEMUDY, DEPDT, DEPDN, DEPDY,
     3  DESDT, DESDN, DESDY, DEUDT, DEUDN, DEUDY,
     4  DPPDN, DPPDT, DPPDY, DPSDN, DPSDT, DPSDY,
     5  DPUDN, DPUDT, DPUDY
C

C
C

