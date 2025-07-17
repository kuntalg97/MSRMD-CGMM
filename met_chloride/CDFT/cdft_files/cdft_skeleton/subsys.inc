!  List of pseudopotentials : https://github.com/cp2k/cp2k/blob/master/data/POTENTIAL
   &SUBSYS
    &CELL
      ABC 20 20 20
      PERIODIC  XYZ
      MULTIPLE_UNIT_CELL  1 1 1
    &END CELL
    &KIND C
      BASIS_SET TZV2P-GTH
      POTENTIAL GTH-BLYP-q4
    &END KIND
    &KIND H
      BASIS_SET TZV2P-GTH
      POTENTIAL GTH-BLYP-q1
    &END KIND
    &KIND O
      BASIS_SET TZV2P-GTH
      POTENTIAL GTH-BLYP-q6
    &END KIND
    &KIND Cl
      BASIS_SET TZV2P-GTH
      POTENTIAL GTH-BLYP-q7
    &END KIND
     &TOPOLOGY
       !COORD_FILE_NAME input.pdb
       !COORD_FILE_FORMAT  PDB
       COORD_FILE_NAME  ${XYZFILE}
       COORD_FILE_FORMAT  XYZ
!       NUMBER_OF_ATOMS  ${NATOMS} ! not used for xyz file 
       MULTIPLE_UNIT_CELL  1 1 1
     &END TOPOLOGY
   &END SUBSYS
   &PRINT
    @IF ( ${JOB_TYPE} == ENERGY_FORCE )
      &FORCES LOW
         FILENAME force
      &END FORCES
    @ENDIF
    @IF ( ${JOB_TYPE} /= ENERGY_FORCE )
      &FORCES OFF
      &END FORCES
    @ENDIF
    &GRID_INFORMATION
    &END
   &END PRINT
