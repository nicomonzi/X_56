$ $Header$
$ ALTER for MSC/NASTRAN 2024.1 SOL 106
$ SOL 106 NONLINEAR STATIC + MODAL ANALYSIS 
$
$ MBDyn (C) is a multibody analysis code.
$ http://www.mbdyn.org
$ 
$ USAGE:
$   In your BDF file:
$     INCLUDE 'MBDyn_NASTRAN_alter_SOL106_2024.nas'
$     CEND
$
$   Enable standard .op2 output in CASE CONTROL:
$     DISPLACEMENT(PLOT) = ALL
$     VECTOR(PLOT) = ALL
$
$   Enable .op2 output with geometry in the BULK DATA:
$     PARAM,OGEOM,YES
$     PARAM,POST,-1
$
$ OUTPUT FILES:
$   mbdyn_modal.mat : Modal matrices (MHH, KHH, LUMPMS) via OUTPUT4
$   <model>.op2     : Geometry (GEOM1), static displacements, and mode shapes (OUG1)
$
$ RUNNING FEMGEN:
$   femgen mymodel              -> reads mymodel.op2, mbdyn_modal.mat -> mymodel.fem
$   femgen mymodel -o output    -> reads mymodel.op2, mbdyn_modal.mat -> output.fem
$
$ NOTE: femgen automatically finds mbdyn_modal.mat if <model>.mat doesn't exist
$
$--------------------------------------------------------------------------
ASSIGN OUTPUT4='mbdyn_modal.mat' STATUS=UNKNOWN UNIT=15
SOL 106


TIME 500

$ Output modal matrices from MODERS subdmap after XREAD
COMPILE      MODERS  
ALTER        124
$ Rename and output modal mass and stiffness
EQUIVX       MIX/MHH/-1 $
LAMX         , ,LAMA/KHH/-1 $
OUTPUT4      MHH,KHH,,//-1/15 $
ENDALTER

$ Output lumped mass from SEDRCVR (continue writing, file stays open until job ends)
COMPILE      SEDRCVR
ALTER        2700
DIAGONAL     MGG/LUMPMS/'COLUMN'/1.  $
OUTPUT4      LUMPMS,,,//-2/15 $
ENDALTER

