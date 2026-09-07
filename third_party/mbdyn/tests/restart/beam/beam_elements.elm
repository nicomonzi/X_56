body: 1, node_id_1 + number_of_nodes - 1,
100.,
null,
diag, 1., 1., 1.;

set: i = 1;
include: "beam_i.elm";

set: i = 2;
include: "beam_i.elm";

set: i = 3;
include: "beam_i.elm";

set: i = 4;
include: "beam_i.elm";

set: i = 5;
include: "beam_i.elm";

set: i = 6;
include: "beam_i.elm";

set: i = 7;
include: "beam_i.elm";

set: i = 8;
include: "beam_i.elm";

set: i = 9;
include: "beam_i.elm";

set: i = 10;
include: "beam_i.elm";

set: i = 1;
include: "body_i.elm";

set: i = 2;
include: "body_i.elm";

set: i = 3;
include: "body_i.elm";

set: i = 4;
include: "body_i.elm";

set: i = 5;
include: "body_i.elm";

set: i = 6;
include: "body_i.elm";

set: i = 7;
include: "body_i.elm";

set: i = 8;
include: "body_i.elm";

set: i = 9;
include: "body_i.elm";

set: i = 10;
include: "body_i.elm";

set: i = 11;
include: "body_i.elm";

set: i = 12;
include: "body_i.elm";

set: i = 13;
include: "body_i.elm";

set: i = 14;
include: "body_i.elm";

set: i = 15;
include: "body_i.elm";

set: i = 16;
include: "body_i.elm";

set: i = 17;
include: "body_i.elm";

set: i = 18;
include: "body_i.elm";

set: i = 19;
include: "body_i.elm";

set: i = 20;
include: "body_i.elm";

set: i = 21;
include: "body_i.elm";

        joint: joint_id_clamp, clamp,
                node_id_1, node, node;

        force: force_id_1, absolute, node_id_1 + number_of_nodes - 1,
                position, reference, ref_node_1 + number_of_nodes - 1, null,
                0., -1., 0.,
                string, "(Time - initial_time) / (final_time - initial_time) * F1";
  gravity: uniform, gx, gy, gz, string, "(Time - initial_time) / (final_time - initial_time)";