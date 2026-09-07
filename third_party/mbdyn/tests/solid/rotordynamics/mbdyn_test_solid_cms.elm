## cond(Mred)=1.7e+03
## cond(Kred)=2.5e+04

joint: elem_id_rotor, modal, node_id_rotor,
        22,
        from file,
        damping from file,
        "mbdyn_test_solid_cms.fem",
        origin node, 654,
        2,
                655, node_id_bearing1, null,
                656, node_id_bearing2, null;
