# X-56A 10 lb modal-shape comparison

This case compares the free-free modes from the 10 lb Nastran SOL 103 model
with the same shapes as reconstructed by the MBDyn modal joint.

## Directory contents

- `nastran/sol103_10lb_f06.bdf`: MSC Nastran SOL 103, 40 modes.
- `nastran/MBDyn_ALTER_103.nas`: writes `mbdyn_modal.mat` for `femgen`.
- `nastran/rbe3s.bdf`: local copy of the established interpolation-grid file.
- `nastran/BULK/`: complete local copy of the released structural model.
- `mbdyn/modal_mode_check.mbd`: activates and excites one flexible modal DOF.
- `compare_mode_shapes.py`: compares the Nastran F06 vectors with the modal
  shapes in the `.fem` file along the signed span.

The structural include set is the released configuration 24611 with
`MAD_L108_FUEL_R11_LT10_NASA.dat` and
`MAD_L108_FUEL_R11_RT10_NASA.dat`. It is the same model used by
`NASTRAN/NASTRAN40/MAIN/nsvibe_test.bdf`.

## Run Nastran on the other machine

Copy only the complete `MODAL_COMPARISON_10LB/nastran` directory. All BDF,
ALTER, RBE3, and structural BULK dependencies are contained inside it.

Then run from `MODAL_COMPARISON_10LB/nastran`:

```bash
chmod +x run_nastran.sh
NASTRAN_CMD=/path/to/nastran ./run_nastran.sh
```

The important result to copy back is `sol103_10lb_f06.f06`. The deck uses
`PRINT,PLOT`, so eigenvectors are written both to the F06 and to the OP2.
It also produces `mbdyn_modal.mat`; run the established `femgen` workflow on
the OP2/MAT pair if the modal basis is regenerated.

## Run the MBDyn check

The current modal joint source is
`NASTRAN/FEMGEN40/mbdyn_modal.fem` (identical to
`bbf/INCLUDE/mbdyn_modal.fem`). In `mbdyn/modal_mode_check.mbd`, set
`MODE_TO_EXCITE` from 7 through 18 and run:

```bash
chmod +x run_mbdyn.sh
./run_mbdyn.sh
```

Only the chosen modal DOF is active and receives a unit generalized force.
Consequently, the physical displacement of every interface node is
proportional to that single `.fem` mode shape. The short 0.02 s run is a
reconstruction check, not a transient-response analysis.

## Generate comparison diagrams

After copying the new F06 back into `nastran/`:

```bash
python3 compare_mode_shapes.py
```

The script writes one English-labelled PNG per mode and
`plots/modal_comparison_summary.csv`. Each shape is normalized by its maximum
translational vector amplitude. The MBDyn sign is aligned to Nastran before
computing the translational MAC, RMS difference, and maximum nodal difference.
The horizontal coordinate is signed elastic-axis arc length, so the swept and
vertical winglet segment is retained.

To compare the actual physical-node reconstruction from a single MBDyn run
instead of reading that shape directly from the `.fem` file:

```bash
python3 compare_mode_shapes.py --modes 7 \
  --mbdyn-mov mbdyn/mode_check.mov --mov-mode 7
```

Before the new printed F06 is available, all direct comparisons between the
existing Nastran modal export and actual MBDyn MOV reconstructions can be
generated with:

```bash
python3 generate_mbdyn_comparisons.py
```

These diagrams are written to `comparison_results/`. The same command also
creates one MBDyn 3D deformed-shape plot per mode. It overlays the light-gray
undeformed line and blue deformed line, uses equal `X`, `Y`, `Z` scales in
inches, and shows the model in a near-isometric view with positive X receding
into the screen. Axis limits follow the geometry, so the Z extent remains
compact. The default displayed maximum modal displacement is 25 inches and
can be changed with `--display-amplitude`. Plot titles include the mode number,
frequency, standard mode acronym, and English mode description.

## Arrange a 3D plot interactively

To arrange one plot manually without saving anything automatically:

```bash
python3 interactive_modal_plot.py --mode 7
```

Drag on the plot to rotate it and drag the legend to place it anywhere. Use
the arrow keys to move the complete graph, the mouse wheel or `+`/`-` to
change zoom, `[`/`]` for line thickness, and `,`/`.` for marker size. Press
`s` only when the desired layout is ready; this writes a timestamped PNG to
`manual_screenshots/`. The Matplotlib toolbar Save button can instead be used
to select the destination and filename.
