        beam3: beam_id_beam_1 + i - 1,
                # node 1
                node_id_1 + 2 * ( i - 1), position, reference, ref_node_1 + 2 * ( i - 1), null,
                orientation, reference, ref_node_1 + 2 * ( i - 1), eye,

                # node 2
                node_id_1 + 2 * ( i - 1 ) + 1, position, reference, ref_node_1 + 2 * ( i - 1) + 1, null,
                orientation, reference, ref_node_1 + 2 * ( i - 1) + 1, eye,

                # node 3,
                node_id_1 + 2 * ( i - 1 )  + 2, position, reference, ref_node_1 + 2 * ( i - 1) + 2, null,
                orientation, reference, ref_node_1 + 2 * ( i - 1) + 2, eye,

                # orientation matrix section I
                reference, ref_gauss_12_1 + i - 1 , eye,
                # reference, ref_node_1 + 2 * ( i - 1) , eye,

                # constitutive law section I
                reference, const_law_id_beam,

                # orientation matrix section II
                reference, ref_gauss_23_1 + i - 1, eye,
                # reference, ref_node_1 + 2 * ( i - 1) +2, eye,

                # constitutive law section II
                reference, const_law_id_beam;
