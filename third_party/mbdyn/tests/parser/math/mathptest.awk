BEGIN {
	if (MBDYN == "") {
		MBDYN = "mbdyn";
	}
}
/^[^#]/ {
	escaped_1 = $1;
	gsub("\"", "\\\"", escaped_1);
	print "echo 'remark: \x22[" $3 "] " escaped_1 " = " $2 "\x22, " $1 ";' | " MBDYN " | grep 'line' | sed -e 's/^[^:]\x5c+: //'";
}
