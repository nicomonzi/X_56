genel: elem_id_m_0 + i,
       spring support,
              node_id_T_0 + i + 1, thermal, algebraic,
              linear viscous generic, m_i * c;

genel: elem_id_kz_0 + i,
       spring,
                node_id_T_0 + i + 1, thermal, algebraic,
                node_id_Tamb, thermal, algebraic,
                linear elastic generic, Aa_i * alpha_side;
