$ ALTER for MSC Nastran 2024.1 SOL 103.
$ Writes the modal matrices required by femgen to mbdyn_modal.mat.
$
$ Standard OP2 geometry/mode-shape output is requested in the parent BDF.

ASSIGN OUTPUT4='mbdyn_modal.mat' STATUS=UNKNOWN UNIT=15

SOL 103
TIME 500

$ Output modal mass and stiffness matrices after XREAD.
COMPILE      XREAD
ALTER        131
EQUIVX       MIX/MHH/-1 $
LAMX         , ,LAMA/KHH/-1 $
OUTPUT4      MHH,KHH,,//-1/15 $
ENDALTER

$ Output the diagonal lumped mass matrix.
COMPILE      SEDRCVR
ALTER        2700
DIAGONAL     MGG/LUMPMS/'COLUMN'/1.  $
OUTPUT4      LUMPMS,,,//-2/15 $
ENDALTER
